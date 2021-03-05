

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


def generateRippsfromTranscripts(transcripts):
	dicRippsPerTranscript = {}
	def findCyclicPermutation(a):
		return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
	def checkMass(nrpseq):
		if len(nrpseq) > 12:
			return 
		if len(nrpseq) < 5:
			return 
		mass = round(sum(nrpseq),2)

		if mass >500 and mass <1800:
			if allSpecMasses[mass] == 1:
				final_nrps_info[mass].append(tuple(findCyclicPermutation((nrpseq))))
		for transcript in transcripts:
			dicRippsPerTranscript[transcripts] = []




def generateNRPSpredict(allContigsBGCdomains,allContigsBGCdomainsScore,k,allSpecMasses):
	allContigsFinalNRP = {}
	import itertools
	for contig in allContigsBGCdomains:
		bgcDomains = allContigsBGCdomains[contig]

		bgc_nrps= {}
		bgc_nrps[0] = {0:[]}
		n = 0
		checked = {}
		# if len(bgcDomains) == 1:
		# 		allContigsFinalNRP[contig] = [[x] for x in  bgcDomains]
		allContigsFinalNRP[contig] = [x[0] for x in itertools.product(itertools.product(*bgcDomains.values()))]


	final_nrps_info = {}
	for possibleMass in [round(p/100.0,2) for p in range(50000, 180000)]:
		final_nrps_info[round(possibleMass,2)] = []

	def checkSingleAddition(orig_nrp):
		orig_mass = round(sum(orig_nrp),2)
		if orig_mass >450 and orig_mass <1800:
			for i in range(len(orig_nrp)):
				temp_full_nrp = list(orig_nrp)

				modified_mass = round(orig_mass + orig_nrp[i],2)
				if modified_mass>500 and modified_mass <1800:
					if allSpecMasses[modified_mass] == 1:
							a= temp_full_nrp[i]		
							temp_full_nrp.insert(i, a)
							modified_nrp = tuple(temp_full_nrp)
							final_nrps_info[modified_mass].append(tuple(modified_nrp))

	def checkSingleDeletion(orig_nrp):
		orig_mass = round(sum(orig_nrp),2)
		if orig_mass >500 and orig_mass <1800:
			for i in range(len(orig_nrp)):
				temp_full_nrp = list(orig_nrp)

				modified_mass = round(orig_mass - orig_nrp[i],2)
				if modified_mass>500 and modified_mass <1800:
					if allSpecMasses[modified_mass] == 1:						
							del temp_full_nrp[i]
							modified_nrp = tuple(temp_full_nrp)
							checkDoubleDeletion(modified_nrp)				
							final_nrps_info[modified_mass].append(tuple(modified_nrp))
	def checkDoubleDeletion(orig_nrp):
		orig_mass = round(sum(orig_nrp),2)
		if orig_mass >500 and orig_mass <1800:
			for i in range(len(nrp)):
				single_modified_mass = round(orig_mass - orig_nrp[i],2)
				if single_modified_mass>500 and single_modified_mass <1800:
					single_modified_nrp = list(orig_nrp)	
					del single_modified_nrp[i]
					for j in range(len(orig_nrp)-1):
						temp_singleDel_nrp = list(single_modified_nrp)
						double_modified_mass = round(single_modified_mass - temp_singleDel_nrp[i],2)
						if double_modified_mass>500 and double_modified_mass<1800:
							if allSpecMasses[double_modified_mass] == 1:						
									del temp_singleDel_nrp[i]
									double_modified_nrp = tuple(temp_singleDel_nrp)					
									final_nrps_info[double_modified_mass].append(tuple(double_modified_nrp))

	from collections import deque
	def findCyclicPermutation(a):
		# for x in sorted([ a[n:] + a[:n] for n in range(len(a)) ]):
		# 	print x
		return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
	def checkMass(nrpseq):
		if len(nrpseq) > 12:
			return 
		if len(nrpseq) < 5:
			return 
		# print nrpseq
		if nrpseq == (129.042, 113.084, 113.084, 99.068, 115.027, 113.084, 113.084):
			print "=============================="
			print "=============================="
			print nrpseq
			print round(sum(nrpseq),2)
		mass = round(sum(nrpseq),2)

		if mass >500 and mass <1800:
			if allSpecMasses[mass] == 1:

				# final_nrps_info[mass].append(tuple(nrpseq))
				final_nrps_info[mass].append(tuple(findCyclicPermutation((nrpseq))))
			

	def mix2Contigs(c1,c2):
		mixedNRP = []
		[checkMass(nrp1 + nrp2) for nrp1, nrp2 in itertools.product(allContigsFinalNRP[c1],allContigsFinalNRP[c2])]
		# for nrp1 in allContigsFinalNRP[c1].values():
		# 	for nrp2 in allContigsFinalNRP[c2].values():
		# 		checkMass(nrp1 + nrp2)

		# return mixedNRP
	def mix3contigs(c1,c2,c3):
		mixedNRP = []
		[checkMass(nrp1 + nrp2 + nrp3) for nrp1, nrp2, nrp3 in itertools.product(allContigsFinalNRP[c1],allContigsFinalNRP[c2],allContigsFinalNRP[c3])]
		# for nrp1 in allContigsFinalNRP[c1].values():
		# 	for nrp2 in allContigsFinalNRP[c2].values():
		# 		for nrp3 in allContigsFinalNRP[c3].values():
		# 			checkMass(nrp1 + nrp2 + nrp3)
					# mixedNRP.append(nrp1 + nrp2 + nrp3)
		# return mixedNRP
	def mix4Contigs(c1,c2,c3,c4):
		mixedNRP = []
		for nrp1 in allContigsFinalNRP[c1].values():
			for nrp2 in allContigsFinalNRP[c2].values():
				for nrp3 in allContigsFinalNRP[c3].values():
					for nrp4 in allContigsFinalNRP[c4].values():
						checkMass(nrp1 + nrp2 +nrp3 + nrp4)
						# mixedNRP.append(nrp1 + nrp2 + nrp3 + nrp4)
		# return mixedNRP
	# allNRPSMIXED = []
	# for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values())):
	# 	print tuple(sum(list(*nrp), () ))
	
	########generate a linear mixture of all ORFs!
	# [checkMass(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a linear mixture!
	orfnames_sorted = sorted(allContigsFinalNRP.keys())[0:min(3,len(allContigsFinalNRP))]
	# for x in [sum(list(*nrp), () ) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in orfnames_sorted]) ) ]:
	# 	print x
	# exit()
	# for z in itertools.product(*[allContigsFinalNRP[orf] for orf in orfnames_sorted]):
	# 	print z
	
	# print sum(list(*nrp), () )
	# exit()
	# print sum(list(tuple(*nrp)), ())
	[checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in orfnames_sorted]) ) ]  #generate a linear mixture!
	########generate a single deletion for each linear mixture of all ORFs
	[checkSingleDeletion(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion
	###### generate a double deletion for each linear mixture of all ORFs
	[checkDoubleDeletion(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion
	######## generate a candidate for each  ORF:
	# [checkMass(nrp) for nrp in [allContigsFinalNRP[contig] for contig in allContigsFinalNRP]] #single contig
	######## generate a candidate for each combination of two ORFs:
	# [checkMass(tuple(list(nrp[0])+list(nrp[1]))) for nrp in itertools.product(itertools.combinations(allContigsFinalNRP.values(),2)) ]
	####### generatea a single amino acid repeat for each linear mixture of all ORFs
	[checkSingleAddition(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion
	
	print "Heeeeeeerekkfjlsjfl"
	# print allContigsFinalNRP
	# print allContigsFinalNRP.values()
	# for nrp in itertools.combinations(allContigsFinalNRP.values(),2):
	# 	print nrp
	# 	print tuple(list(nrp[0])+list(nrp[1]))
	# print allContigsFinalNRP
	# for pair in [itertools.product(allContigsFinalNRP[pair[0]],allContigsFinalNRP[pair[1]]) for pair in itertools.combinations(allContigsFinalNRP.keys(),2) ]:
	# 	print pair
	# for pair in itertools.combinations(allContigsFinalNRP.keys(),2):
	# 	for npr in itertools.product(allContigsFinalNRP[pair[0]],allContigsFinalNRP[pair[1]]):
	# 		checkMass(tuple(list(npr[0])+list(npr[1])) )

	# [itertools.product(allContigsFinalNRP[pair[0]],allContigsFinalNRP[pair[1]]) for pair in itertools.combinations(allContigsFinalNRP.keys(),2) ] ]
	# print [nrp[0] for nrp in [ for pair in itertools.comb inations(allContigsFinalNRP.keys(),2) ]

	# [checkMass(nrp[0]) for nrp in itertools.combinations(allContigsFinalNRP.values(),2) ]
	# if len(allContigsFinalNRP) == 2:
	# 	contig1,contig2 = allContigsFinalNRP.keys()
	# 	mix2Contigs(contig1,contig2)
	# 	# allNRPSMIXED += mix2Contigs(contig2,contig1)
	# 	# for contig1 in allContigsFinalNRP:
	# 	# 	for contig2 in allContigsFinalNRP:
	# 	# 		if contig1 == contig2:
	# 	# 			continue
	# 	# 		allNRPSMIXED += mix2Contigs(contig1,contig2)
	# if len(allContigsFinalNRP) == 3:
	# 	contig1,contig2,contig3 = allContigsFinalNRP.keys()
	# 	# contig1 = sorted(allContigsFinalNRP.keys())[0]
	# 	# contig2 = sorted(allContigsFinalNRP.keys())[1]
	# 	# contig3 = sorted(allContigsFinalNRP.keys())[2]
	# 	mix3contigs(contig1,contig2,contig3)
	
	# if len(allContigsFinalNRP) == 4:
	# 	contig1 = sorted(allContigsFinalNRP.keys())[0]
	# 	contig2 = sorted(allContigsFinalNRP.keys())[1]
	# 	contig3 = sorted(allContigsFinalNRP.keys())[2]
	# 	contig4 = sorted(allContigsFinalNRP.keys())[3]
		# allNRPSMIXED += mix4Contigs(contig1,contig2,contig3,contig4)
	#print len(allNRPSMIXED)

		
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

	
	# for nrps in allNRPSMIXED:
	# 	# allkmers = generateNRPS3mers(nrps,k)
	# 	# allkmers = set()
	# 	mass = round(sum(nrps),2)
	# 	if mass >500 and mass <1800:
	# 		if allSpecMasses[mass] == 1:
	# 			final_nrps_info[mass].append(tuple(nrps))
		
	# temp_final_nrps_info = final_nrps_info.copy()

	# #create nrps with one Adomain deleted
	# for mass in temp_final_nrps_info:
	# 	if len(temp_final_nrps_info[mass])==0:
	# 		continue
	# 	for nrps in temp_final_nrps_info[mass]:
	# 		for j in range(len(nrps)-1):
	# 			mass_new_nrps = round(mass - nrps[j],2)
	# 			if mass_new_nrps > 500 and mass_new_nrps < 1800:
	# 				if allSpecMasses[mass_new_nrps] == 1:
	# 					new_nrps = nrps[:j] + nrps[j+1:]
	# 				else:
	# 					continue
	# 			else:
	# 				continue
	# 			# if tuple(new_nrps) in final_nrps_info:
	# 			# 	continue
	# 	#		allkmers = generateNRPS3mers(new_nrps,k)
	# 			# allkmers = set()
	# 			# final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	# 			# if allSpecMasses[mass_new_nrps] == 1:
	# 				# final_nrps_info[round(sum(new_nrps),2)].update({tuple(new_nrps):round(sum(new_nrps),5)})
	# 			final_nrps_info[round(sum(new_nrps),2)].append(tuple(new_nrps))

	


	# final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]
	# 	new_nrps = nrps[:len(nrps)-1]
	# 	allkmers = generateNRPS3mers(new_nrps,k)
	# 	final_nrps_info[tuple(new_nrps)]= [round(sum(new_nrps),5),set(allkmers)]


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

	# candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
	# candidateNRPSs_sorted = candidateNRPSs.items()

	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False
	foundOne = False
	for nrpsscore in candidateNRPSs:
		nrps = nrpsscore
		# score = nrpsscore[1]
		# score = 3
		# if score>2:
		graphFile = open(output+"_"+str(precursorMass)+"_candidateNRPs.GRAPHS.txt","a")
		linesToAdd = createPNPGraphLines(nrps)
		foundOne = True
		for line in linesToAdd:
			graphFile.write(line+"\n")
	return foundOne
		# arguments_to_print = [",".join([str(x) for x in nrps]),precursorMass, retention, charge, score]
		# reconstructions_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")


def output_nrpsgraph(candidateNRPSs, output):
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

	# candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
	# candidateNRPSs_sorted = candidateNRPSs.items()

	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False
	foundOne = False
	for nrpsscore in candidateNRPSs:
		nrps = nrpsscore
		# score = nrpsscore[1]
		# score = 3
		# if score>2:
		graphFile = open(output+"_candidateNRPs.GRAPHS.txt","a")
		linesToAdd = createPNPGraphLines(nrps)
		foundOne = True
		for line in linesToAdd:
			graphFile.write(line+"\n")
	return foundOne

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
