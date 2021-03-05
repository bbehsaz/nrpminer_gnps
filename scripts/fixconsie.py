import sys, string
n=0
#fixedfile = open("fixed_"+sys.argv[1],"w")
inputfile = open(sys.argv[1],"r")
k = int(sys.argv[2])
def fixConsivePrintScoreOutput(intpufile):
	#inputfile = open(inputfilename)
	line2print=""
	while True:
		line = inputfile.readline()
		# print line
		if not line:
			break
		if (n % k) == 0:
			line2print = line.strip() + "\t"
			for i in range(k-1):
				line=inputfile.readline()
				line2print+=line.strip().replace(" ","_")+"\t"
			print line2print
fixConsivePrintScoreOutput(inputfile)