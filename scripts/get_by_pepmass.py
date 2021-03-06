import sys, string

import sys, string

def readStream2(peaksFile):

	protonMass = 1.00728
	peaks = {}
	peaksnIntensity = {}
	pepMasses={}
	retention=-1 	
	charges= {}
	n= -1
	retentions = {}
	specLines = {}
	fileLines = {}
	#Reading the peaks and intensities.
	while(True):
		line = peaksFile.readline()
		if not line:
			break
		originalLine = line
		line = line.strip()
		if line == "BEGIN IONS":
			specLines = [originalLine]
			peptide= ""			
			while(True):
				line = peaksFile.readline()
				if not line[0].isdigit():
					specLines.append(line)
					if line[0:6]=="CHARGE":
						charge = int(''.join(c for c in line.split("=")[1] if c.isdigit()))						
					if line[0:6]=="PEPMAS":
						pepMass=float(line.strip().split()[0][8:])
					if line[0:5]=="DBID=":
						line= line.strip()
						if line[5:] !="":
							antimartin = line[15:]
					if "TITLE=" in line:
						peptide=line.strip().split("=")[1]
					if "RTINSECONDS" in line:
						retention = line.split("=")[1].strip()
				else:
					break
			n +=1
			if peptide == "":
				peptide=str(n)
			peptide = peptide + "_"+str(n) 
			retentions[peptide] = retention
			peaksnIntensity[peptide]={}
			charges[peptide]=charge
			pepMasses[peptide]= pepMass
			while line!="END IONS":
				specLines.append(line.strip()+"\n")
				peakLine = line.strip().split()
				peakMass = round(float(peakLine[0]),3)
				intensity = float(peakLine[1])
				peaksnIntensity[peptide][peakMass]=intensity
				if charge == 2:
					charge1PeakMass = (2*peakMass) - protonMass
					peaksnIntensity[peptide][charge1PeakMass]=intensity
				line = peaksFile.readline().strip()
			if line == "END IONS":
				specLines.append(line+"\n")
				fileLines[peptide] = specLines
				continue
		
	return peaksnIntensity,pepMasses,charges,retentions,fileLines

streamFile = sys.stdin

peaksnIntensity,pepMasses,charges,retentions,fileLines= readStream2(streamFile)

masses = set()

massfile = open(sys.argv[1])
while True: #first file is a file containing all exact masses per each line
	line = massfile.readline()
	if not line:
		break
	#print "ELEAVE ME THE FUCK"
	#print line
	masses.add(float(line.strip()))
	#print "wHATTHEFUCK"
	#print line.strip()

for peptide in peaksnIntensity:
	#if pepMasses[peptide] == float(sys.argv[1]):
	if pepMasses[peptide] in masses:
		outputfile = open(sys.argv[2]+str(pepMasses[peptide])+".mgf",'a') #second argument output dir
		for line in fileLines[peptide]:
			outputfile.write(line)


