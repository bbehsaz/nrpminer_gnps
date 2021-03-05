import sys, string

#THIS is right now only for one peptide


def getBuildingBlocksForNRPSprediction(path2bbFile): 
#The function gets the building block files and return the required datastructures 
	digits = 3
	aminoMAsses = {}
	standardMasses = set()
	mass2name = {}
	aminFullName2Mass = {}
	with open(path2bbFile) as bbFile:
		for line in bbFile:
				linesplit = line.strip().split()
				aminoMAsses[linesplit[0]] = round(float(linesplit[3]),digits)
				aminFullName2Mass[linesplit[1].lower()] = round(float(linesplit[3]),digits)
				standardMasses.add(round(float(linesplit[3]),digits))
				mass2name[round(float(linesplit[3]),digits)] = linesplit[0]
	return standardMasses,aminFullName2Mass

def gbkFile(gbkFileAddress,topAA): #we do not filter at this point on any point. JUST READ the codes file.
	# aaMassesFile = open("./configs/aminoacidMasses.txt","r")
	# aaMassesFile.readline()
	# aaMasses = {}
	# while(True):
	# 	line = aaMassesFile.readline()
	# 	if not line:
	# 		break
	# 	aaMasses[line.split()[0]] = float(line.split()[5])
	standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction('./configs/nrp_hosein_mass.txt')
	antiSmashPredicts ={}
	antiSmashPredictsScore = {}
	
	gbkFile = open(gbkFileAddress,"r")
	line =gbkFile.readline().strip()
	new_contig = line.split()[0].split("_A")[0] #these are orfs not contigs
	while(True):
		if not line:
			break
		lineSplit = line.split()
		predictions = line.split()[2]
		aminoAcids = predictions.split(";")
		contig = new_contig
		pos = 0
		antiSmashPredicts[contig] = {}
		antiSmashPredictsScore[contig] = {}
		while new_contig ==contig:
			antiSmashPredicts[contig][pos] = []
			antiSmashPredictsScore[contig][pos] = []
			
			for i in range(15):
				if len(set(antiSmashPredictsScore[contig][pos]))>topAA-1:
					continue
				aminoScore = aminoAcids[i]
				aminoName = aminoScore.split("(")[0]
				score = float(aminoScore.split("(")[1].split(")")[0])
				if aminoName in aminFullName2Mass:
					if score<60:
						continue
					antiSmashPredicts[contig][pos].append(aminFullName2Mass[aminoName])
					antiSmashPredictsScore[contig][pos].append(score)
			pos += 1
			line =gbkFile.readline().strip()
			if not line: 
				break
			lineSplit = line.split()
			new_contig =  lineSplit[0].split("_A")[0]
			predictions = line.split()[2]
			aminoAcids = predictions.split(";")
	finalAntiSmashPredicts = {}
	finalAntiSmashPredictsScore = {}
	for contig in antiSmashPredicts:
		if len(antiSmashPredicts[contig]) <0:
			continue
		finalAntiSmashPredicts[contig] = antiSmashPredicts[contig]
		finalAntiSmashPredictsScore[contig] = antiSmashPredictsScore[contig]
		# print contig
		# print antiSmashPredicts[contig]
	return finalAntiSmashPredicts,finalAntiSmashPredictsScore


def generateNRPS3mers(nrps,k):
	def creatkmername(kmertuple):
		return "-".join([str(x) for x in kmertuple])
	nrpskmers = {}
	for i in range(len(nrps)):
		if i+k<len(nrps)+1:
			kmer = tuple(nrps[i:i+k])
		else:
			kmer = tuple(nrps[i:]+nrps[:k-(len(nrps)-i)])		
		nrpskmers[creatkmername(kmer)] = 1
		reverseKmer = kmer[::-1]

		nrpskmers[creatkmername(reverseKmer)] = 1
	return nrpskmers


def generateNRPSpredict(allContigsBGCdomains,allContigsBGCdomainsScore,k):
	allContigsFinalNRP = {}
	for contig in allContigsBGCdomains:
		bgcDomains = allContigsBGCdomains[contig]
		bgc_nrps= {}
		bgc_nrps[0] = {0:[]}
		n = 0
		checked = {}
		for i in range(1,len(bgcDomains)+1):
			bgc_nrps[i] = {}
			for sequenceNum in bgc_nrps[i-1]:
				sequence = bgc_nrps[i-1][sequenceNum][:]
				for amino in bgcDomains[i-1]:
					extended_sequence = sequence + [amino]
					if tuple(extended_sequence) not in checked:
						n+=1 
						bgc_nrps[i][n]= extended_sequence[:]
						checked[tuple(extended_sequence)]=1
		allContigsFinalNRP[contig] = bgc_nrps[len(bgcDomains)].copy()
	for contig in allContigsFinalNRP:
		print contig
		print allContigsFinalNRP[contig]
	def mixContigs(c1,c2):
		mixedNRP = []
		for nrp1 in allContigsFinalNRP[c1].values():
			for nrp2 in allContigsFinalNRP[c2].values():
				mixedNRP.append(nrp1 + nrp2)
		return mixedNRP
	allNRPSMIXED = []
	for contig1 in allContigsFinalNRP:
		for contig2 in allContigsFinalNRP:
			if contig1 == contig2:
				continue
			allNRPSMIXED += mixContigs(contig1,contig2)
	final_nrps_info = {}
		
		# #create nrps with one Adomain deleted
		# for nrpsNum in bgc_nrps[len(bgcDomains)]:
		# 	nrps = bgc_nrps[len(bgcDomains)][nrpsNum]
		# 	for j in range(len(bgcDomains)-1):
		# 		new_nrps = nrps[:j] + nrps[j+1:]
		# 		allkmers = generateNRPS3mers(new_nrps,k)
		# 		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		# 	new_nrps = nrps[:len(bgcDomains)-1]
		# 	allkmers = generateNRPS3mers(new_nrps,k)
		# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		# temp_final_nrps_info = final_nrps_info.copy()
		# #create nrps with two Adomain deleted	
		# for nrps in temp_final_nrps_info:
		# 	for j in range(len(nrps)-1):
		# 		new_nrps = nrps[:j] + nrps[j+1:]
		# 		if tuple(new_nrps) in final_nrps_info:
		# 			continue
		# 		allkmers = generateNRPS3mers(new_nrps,k)
		# 		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		# 	new_nrps = nrps[:len(nrps)-1]
		# 	allkmers = generateNRPS3mers(new_nrps,k)
		# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		
	for nrps in allNRPSMIXED:
		allkmers = generateNRPS3mers(nrps,k)
		final_nrps_info[tuple(nrps)]= [round(sum(nrps),5),set(allkmers)]
	temp_final_nrps_info = final_nrps_info.copy()

	# final_nrps_info = {} #right now it only includes one deletion!
	#create nrps with one Adomain deleted	
	for nrps in temp_final_nrps_info:
		for j in range(len(nrps)-1):
			new_nrps = nrps[:j] + nrps[j+1:]
			if tuple(new_nrps) in final_nrps_info:
				continue
			allkmers = generateNRPS3mers(new_nrps,k)
			final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		new_nrps = nrps[:len(nrps)-1]
		allkmers = generateNRPS3mers(new_nrps,k)
		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]


	temp_final_nrps_info = final_nrps_info.copy()
	# final_nrps_info = {} #right now it only includes duplicateds!
	#create nrps with one Adomain duplicated	
	for nrps in temp_final_nrps_info:
		for j in range(len(nrps)-1):
			new_nrps = nrps[:j+1] + tuple([nrps[j]]) + nrps[j+1:]
			if tuple(new_nrps) in final_nrps_info:
				continue
			allkmers = generateNRPS3mers(new_nrps,k)
			final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		new_nrps = nrps[:len(nrps)-1]
		allkmers = generateNRPS3mers(new_nrps,k)
		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]


	temp_final_nrps_info = final_nrps_info.copy()
	# final_nrps_info = {} #right now it only includes duplicateds!
	#create nrps with two Adomain duplicated	
	for nrps in temp_final_nrps_info:
		for j in range(len(nrps)-1):
			new_nrps = nrps[:j+1] + tuple([nrps[j]]) + nrps[j+1:]
			if tuple(new_nrps) in final_nrps_info:
				continue
			allkmers = generateNRPS3mers(new_nrps,k)
			final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
		new_nrps = nrps[:len(nrps)-1]
		allkmers = generateNRPS3mers(new_nrps,k)
		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]

	# temp_final_nrps_info = final_nrps_info.copy()
	# #create nrps with two A domain deleted	
	# for nrps in temp_final_nrps_info:
	# 	for j in range(len(nrps)-1):
	# 		new_nrps = nrps[:j] + nrps[j+1:]
	# 		if tuple(new_nrps) in final_nrps_info:
	# 			continue
	# 		allkmers = generateNRPS3mers(new_nrps,k)
	# 		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	# 	new_nrps = nrps[:len(nrps)-2]
	# 	allkmers = generateNRPS3mers(new_nrps,k)
	# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	return final_nrps_info 



def output_spectrum_nrps(candidateNRPSs, precursorMass, retention, charge, peptide, reconstructions_file):
	#This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
	from operator import itemgetter
	
	candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False

	for nrpsscore in candidateNRPSs_sorted:
		nrps = nrpsscore[0]
		score = nrpsscore[1]
		arguments_to_print = [",".join([str(x) for x in nrps]),precursorMass, retention, charge, score]
		reconstructions_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")

		# scores.add(score)
	# 	if len(scores)>4:
	# 		break
	# 	elif score >12:
	# 		peptidesUptoCorrectScore.append(reconstruction)
	# 		scores.add(score)
	# 	if score>11:
	# 		reconstructions_file.write("{}\t{}\t{}\t{}\t{}\n".format(",".join([str(x) for x in reconstruction]),precursorMass, retention, charge, score))
	# if len(scores) == 0:
	# 	maxScore = 0
	# else:
	# 	maxScore = max(scores)

	# arguments_to_print = [(kmerSize,kmerThreshold), precursorMass, retention, charge, len(peptidesUptoCorrectScore), len(scores), maxScore]
	# benchmark_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")

	# if maxScore>12:
	# 	return 1
	# else:
	# 	return 0






# bgcDomains,bgcDomainsScore= readGBK('/Users/bahar/workspace/antibiotic-sequencing/code/high_analysis/nrp_quest/cyclonovo/fake_surugamide_codes.gbk')
# generateNRPSpredict(bgcDomains,bgcDomainsScore)
# print generateNRPS3mers(bgcDomains,bgcDomainsScore,3)
#print readGBK(sys.argv[1])[0]
