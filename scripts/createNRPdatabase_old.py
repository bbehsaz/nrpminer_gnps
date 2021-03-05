

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
	# for contig in allContigsFinalNRP:
	# 	print contig
	# 	print allContigsFinalNRP[contig]
	def mix2Contigs(c1,c2):
		mixedNRP = []
		for nrp1 in allContigsFinalNRP[c1].values():
			for nrp2 in allContigsFinalNRP[c2].values():
				mixedNRP.append(nrp1 + nrp2)
		return mixedNRP
	def mix3contigs(c1,c2,c3):
		mixedNRP = []
		for nrp1 in allContigsFinalNRP[c1].values():
			for nrp2 in allContigsFinalNRP[c2].values():
				for nrp3 in allContigsFinalNRP[c3].values():
					mixedNRP.append(nrp1 + nrp2 + nrp3)
		return mixedNRP
	def mix4Contigs(c1,c2,c3,c4):
		mixedNRP = []
		for nrp1 in allContigsFinalNRP[c1].values():
			for nrp2 in allContigsFinalNRP[c2].values():
				for nrp3 in allContigsFinalNRP[c3].values():
					for nrp4 in allContigsFinalNRP[c4].values():
						mixedNRP.append(nrp1 + nrp2 + nrp3 + nrp4)
		return mixedNRP
	allNRPSMIXED = []
	#Shuffling two contigs 
	# for contig1 in allContigsFinalNRP:
	# 	for contig2 in allContigsFinalNRP:
	# 		if contig1 == contig2:
	# 			continue
	# 		allNRPSMIXED += mixContigs(contig1,contig2)
	#linear add
	# for i in range(len(allContigsFinalNRP)):
	print "HEEEERE"
	print allContigsFinalNRP.keys()
	if len(allContigsFinalNRP) == 1:
		for nrp in allContigsFinalNRP[allContigsFinalNRP.keys()[0]]:
			allNRPSMIXED+= nrp
	if len(allContigsFinalNRP) == 2:
		contig1 = sorted(allContigsFinalNRP.keys())[0]
		contig2 = sorted(allContigsFinalNRP.keys())[1]

		allNRPSMIXED += mix2Contigs(contig1,contig2)
		# for contig1 in allContigsFinalNRP:
		# 	for contig2 in allContigsFinalNRP:
		# 		if contig1 == contig2:
		# 			continue
		# 		allNRPSMIXED += mix2Contigs(contig1,contig2)
	if len(allContigsFinalNRP) == 3:
		contig1 = sorted(allContigsFinalNRP.keys())[0]
		contig2 = sorted(allContigsFinalNRP.keys())[1]
		contig3 = sorted(allContigsFinalNRP.keys())[2]
		allNRPSMIXED += mix3contigs(contig1,contig2,contig3)
	
	# if len(allContigsFinalNRP) == 4:
	# 	contig1 = sorted(allContigsFinalNRP.keys())[0]
	# 	contig2 = sorted(allContigsFinalNRP.keys())[1]
	# 	contig3 = sorted(allContigsFinalNRP.keys())[2]
	# 	contig4 = sorted(allContigsFinalNRP.keys())[3]
		# allNRPSMIXED += mix4Contigs(contig1,contig2,contig3,contig4)
	print len(allNRPSMIXED)

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
		# allkmers = generateNRPS3mers(nrps,k)
		allkmers = set()
		final_nrps_info[tuple(nrps)]= [round(sum(nrps),5),set(allkmers)]
	
	temp_final_nrps_info = final_nrps_info.copy()

	final_nrps_info = {} #right now it only includes one deletion!
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


	# temp_final_nrps_info = final_nrps_info.copy()
	# # final_nrps_info = {} #right now it only includes duplicateds!
	# #create nrps with one Adomain duplicated	
	# for nrps in temp_final_nrps_info:
	# 	for j in range(len(nrps)-1):
	# 		new_nrps = nrps[:j+1] + tuple([nrps[j]]) + nrps[j+1:]
	# 		if tuple(new_nrps) in final_nrps_info:
	# 		[s ]	continue
	# 		allkmers = generateNRPS3mers(new_nrps,k)
	# 		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	# 	new_nrps = nrps[:len(nrps)-1]
	# 	allkmers = generateNRPS3mers(new_nrps,k)
	# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]


	# temp_final_nrps_info = final_nrps_info.copy()
	# # final_nrps_info = {} #right now it only includes duplicateds!
	# #create nrps with two Adomain duplicated	
	# for nrps in temp_final_nrps_info:
	# 	for j in range(len(nrps)-1):
	# 		new_nrps = nrps[:j+1] + tuple([nrps[j]]) + nrps[j+1:]
	# 		if tuple(new_nrps) in final_nrps_info:
	# 			continue
	# 		allkmers = generateNRPS3mers(new_nrps,k)
	# 		final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	# 	new_nrps = nrps[:len(nrps)-1]
	# 	allkmers = generateNRPS3mers(new_nrps,k)
	# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]

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


def output_spectrum_nrpsgraph(candidateNRPSs, precursorMass, retention, charge, peptide, output):
	#This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
	from operator import itemgetter
	def createPNPGraphLines(pnp):
		num = len(pnp)
		lines= []
		newline = "number of components : "+ str(num)
		lines.append(newline)
		for i in range(num):
			aa = pnp[i]
			newline = str(i) + " CXHXNX " + str(aa)
			lines.append(newline)
		newline = "number of bonds : " + str(num)
		lines.append(newline)
		for i in range(num-1):
			newline = str(i) + " -NC> " + str(i+1)
			lines.append(newline)
		newline = str(num-1) + " -NC> " + str(0)
		lines.append(newline)
		return lines

	candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False
	foundOne = False
	for nrpsscore in candidateNRPSs_sorted:
		nrps = nrpsscore[0]
		score = nrpsscore[1]
		if score>2:
			graphFile = open(output+"_"+str(precursorMass)+"_candidateNRPs.GRAPHS.txt","a")
			linesToAdd = createPNPGraphLines(nrps)
			foundOne = True
			for line in linesToAdd:
				graphFile.write(line+"\n")
	return foundOne
		# arguments_to_print = [",".join([str(x) for x in nrps]),precursorMass, retention, charge, score]
		# reconstructions_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")



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


def createCompoundGraph(pnp,outputGraphFile):
	from operator import itemgetter
	
	candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False

	for nrpsscore in candidateNRPSs_sorted:
		nrps = nrpsscore[0]
		score = nrpsscore[1]

		# arguments_to_print = [",".join([str(x) for x in nrps]),precursorMass, retention, charge, score]
		# reconstructions_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")




# bgcDomains,bgcDomainsScore= readGBK('/Users/bahar/workspace/antibiotic-sequencing/code/high_analysis/nrp_quest/cyclonovo/fake_surugamide_codes.gbk')
# generateNRPSpredict(bgcDomains,bgcDomainsScore)
# print generateNRPS3mers(bgcDomains,bgcDomainsScore,3)
#print readGBK(sys.argv[1])[0]
