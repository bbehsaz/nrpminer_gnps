

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


def generateNRPSpredict(allContigsBGCdomains,allContigsBGCdomainsScore,k,allSpecMasses,minAdomainscore):
    allContigsFinalNRP = {}
    import itertools
    import math

    for orf in allContigsBGCdomains:
        print "----"
        print orf
        bgcDomains = allContigsBGCdomains[orf]
        print bgcDomains
        allFinalScores = {}
        allFinalRanks = {}
        acceptableBGCDomains = {}
        for pos in bgcDomains.keys():
            acceptableBGCDomains[pos] =[]
            allFinalScores[pos] =  {}
            allFinalRanks[pos] = {}
            sortedScores = sorted(set(allContigsBGCdomainsScore[orf][pos]),reverse=True)
            for i in range(len(bgcDomains[pos])):
                if bgcDomains[pos][i] in allFinalScores[pos]:
                    continue

                if math.ceil(allContigsBGCdomainsScore[orf][pos][i])> 0:
                    acceptableBGCDomains[pos].append(bgcDomains[pos][i])
                    allFinalScores[pos][bgcDomains[pos][i]] = int(allContigsBGCdomainsScore[orf][pos][i])
                    allFinalRanks[pos][bgcDomains[pos][i]] = sortedScores.index(int(allContigsBGCdomainsScore[orf][pos][i]))  
    
        bgc_nrps= {}
        bgc_nrps[0] = {0:[]}
        n = 0
        checked = {}

        # allnrps = [x[0] for x in itertools.product(itertools.product(*bgcDomains.values()))]
        allnrps = [x[0] for x in itertools.product(itertools.product(*acceptableBGCDomains.values()))]
        print "allContigsFinalNRP"
        if len(bgcDomains)>0:
            # print "HEEEEREEE"
            # print "bgcDomains"
            # print bgcDomains
            # print "acceptableBGCDomains"
            # print acceptableBGCDomains
            # print "allFinalScores"
            # print allFinalScores
            # print "allFinalRanks"
            # print allFinalRanks
            # allContigsFinalNRP[orf] = [nrp for nrp in allnrps if ([allFinalRanks[pos][nrp[pos]] for pos in range(len(nrp)) ].count(2) ) < math.ceil(len(nrp)*0.3) ]
            allContigsFinalNRP[orf] = [nrp for nrp in allnrps if ([allFinalRanks[pos][nrp[pos]] for pos in range(len(nrp)) ].count(1) ) < 2 ]
        # print allContigsFinalNRP[orf]

    # exit()
    final_nrps_info = {}
    for possibleMass in [round(p/100.0,2) for p in range(50000, 180000)]:
        final_nrps_info[round(possibleMass,2)] = []

    from collections import deque
    def findCyclicPermutation(a):
        # for x in sorted([ a[n:] + a[:n] for n in range(len(a)) ]):
        #   print x
        return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
            


    #########################   generate a linear mixture of all ORFs!
    ####### We first we find all sequential ORFs (i.e. consecutive ID numbers) then mix them in linear fashion. TO DO: Shuffle them around.
    ###sort the ORFs based on their integer numbers to find the sequential ORFs
    import re
    def atoi(text):
        return int(text) if text.isdigit() else text
    def natural_keys(text):
        return [ atoi(c) for c in re.split('(\d+)', text) ]
    orfnames_sorted = allContigsFinalNRP.keys()[:]
    orfnames_sorted.sort(key=natural_keys) #sorted list of ORF names
    
    #### Group them together based on consecutive numbers
    consectuiveORFgroups = {}
    old_orf = -100
    group_pointer = 0
    for orf in orfnames_sorted:
        current_orf_str = natural_keys(orf)[3]
        if not str(current_orf_str).isdigit():
            print "The name of ORF in the nrpspredictor2 output is not consistent with expected output format!"
            return False
        else:
            current_orf = int(current_orf_str)
        if abs(current_orf - old_orf) <12:
            consectuiveORFgroups[group_pointer].append(orf)
        else:
            group_pointer +=1 
            consectuiveORFgroups[group_pointer] = [orf]
        old_orf = current_orf
    return consectuiveORFgroups, allContigsFinalNRP
    ###### Generate a linear mixture of each consecutive group ... for groups with one member only only nrps generated from that ORF is explored


def findNRPPredicts(consectuiveORFgroups,allSpecMasses,allContigsFinalNRP):
    import itertools
    final_nrps_info = {}
    for possibleMass in [round(p/100.0,2) for p in range(50000, 180000)]:
        final_nrps_info[round(possibleMass,2)] = []
    from collections import deque
    def findCyclicPermutation(a):

        return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
    def checkMass(nrpseq):
        # if len(nrpseq) != 8:
        #     return
        if len(nrpseq) > 14:
            return 
        if len(nrpseq) < 2:
            return 
        found = False
        mass = round(sum(nrpseq),2)
        if mass >500 and mass <1800:
            if allSpecMasses[mass] == 1:
                if found:
                    print "We found it here WWWWEWEWEFDSFKJASLDKJASLKDFJLKJSDF"
                    print tuple(findCyclicPermutation((nrpseq)))
                final_nrps_info[mass].append(tuple(findCyclicPermutation((nrpseq))))

    import itertools

    for group in consectuiveORFgroups:
        print "-=============================="
        print group
        print consectuiveORFgroups[group]
        all_orfs_linear_order = consectuiveORFgroups[group][:]
        print all_orfs_linear_order
        # if all_orfs_linear_order != ['ctg1_orf00658', 'ctg1_orf00659', 'ctg1_orf00660', 'ctg1_orf00664', 'ctg1_orf00665']:
        #     continue
        print "#of actual nrps: {}".format(len([nrp for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]))

        # print [allContigsFinalNRP[orf] for orf in all_orfs_linear_order]

        # print [x for x in itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) ]

        # print [nrp for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]
        # print [(list(*nrp), () ) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]
        [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]


    #   # one_orf_delet_list = [] 
    #   # print [(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]
    #   # map(set, itertools.combinations(all_orfs_linear_order, 2))

    #   ######## create all pairs of two orfs
    #   # for orf_pair in map(list, itertools.combinations(all_orfs_linear_order, 2)):
    #   #   print sorted((orf_pair))
    #   #   [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in sorted(list(orf_pair)) ] ) )]
    #   # print x[0:3]
    #   # y=[(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in all_orfs_linear_order ] ) )]
    #   # print y[0:3]
    #   ###### straight linear version
      # product = 1
      # for x in [len(allContigsFinalNRP[orf]) for orf in all_orfs_linear_order ]:
      #     product *= x

      # print "#of actual nrps: {}".format(x)
# 
        ################################ 
        # [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in one_orf_delet_list ] ) )]
        # for one_del_orf_seq in one_orf_delet_list:
        #   [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in one_del_orf_seq] ) )]
        # [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in two_orf_delet_list[0] ] ) ) ]
        # for two_del_orf_seq in two_orf_delet_list:
        #   # print [sum(list(*nrp), () ) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in two_del_orf_seq] ) )]
        #   [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in two_del_orf_seq] ) ) ]

        # [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in consectuiveORFgroups[group] ] ) ) ]  #generate a linear mixture!
    # [checkMass(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in orfnames_sorted]) ) ]  #top 3!!
    ########generate a single deletion for each linear mixture of all ORFs
    # for group in consectuiveORFgroups:
    #   [checkSingleDeletion(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in consectuiveORFgroups[group] ] ) )]
    # [checkSingleDeletion(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion
    # ###### generate a double deletion for each linear mixture of all ORFs
    # for group in consectuiveORFgroups:
    #   [checkDoubleDeletion(sum(list(*nrp), () )) for nrp in itertools.product(itertools.product(*[allContigsFinalNRP[orf] for orf in consectuiveORFgroups[group] ] ) )]
    # [checkDoubleDeletion(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion
    # ######## generate a candidate for each  ORF:
    # # [checkMass(nrp) for nrp in [allContigsFinalNRP[contig] for contig in allContigsFinalNRP]] #single contig
    # ######## generate a candidate for each combination of two ORFs:
    # # [checkMass(tuple(list(nrp[0])+list(nrp[1]))) for nrp in itertools.product(itertools.combinations(allContigsFinalNRP.values(),2)) ]
    # ####### generatea a single amino acid repeat for each linear mixture of all ORFs
    # [checkSingleAddition(tuple(sum(list(*nrp), () ))) for nrp in itertools.product(itertools.product(*allContigsFinalNRP.values()))] #generate a single deletion

    print "Heeeeeeerekkfjlsjfl"
    return final_nrps_info 



def output_spectrum_cyclic_nrpsgraph(candidateNRPSs, precursorMass, retention, charge, peptide, output):
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

    scores = set()
    
    peptidesUptoCorrectScore = []
    correctFound = False
    foundOne = False
    for nrpsscore in candidateNRPSs:
        nrps = nrpsscore
        graphFile = open(output+"_"+str(precursorMass)+"_candidateNRPs.GRAPHS.txt","a")
        linesToAdd = createPNPGraphLines(nrps)
        foundOne = True
        for line in linesToAdd:
            graphFile.write(line+"\n")
    return foundOne


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
def output_nrpsgraph_bcyclic(candidateNRPSs, output):
    #This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
    from operator import itemgetter
    protonMass = 1.00728 #This is simply mass of hydrogen --- do not remove!
    def createPNPGraphLines(pnp):
        num = len(pnp)
        lines= []
        newline = "number of components : "+ str(num)
        lines.append(newline)
        for i in range(num):
            aa = pnp[i]
            if i == 1:
                aa -= protonMass
            newline = str(i) + " CXHXNX " + str(aa)

            lines.append(newline)
        newline = "number of bonds : " + str(num)
        lines.append(newline)
        newline = str(0) + " -NC> " + str(1)
        lines.append(newline)
        for i in range(1,num-1):
            newline = str(i) + " -NC> " + str(i+1)
            lines.append(newline)
        newline = str(num-1) + " -NC> " + str(1)
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



def createCompoundGraph(pnp,outputGraphFile):
    from operator import itemgetter
    
    candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)
    scores = set()
    
    peptidesUptoCorrectScore = []
    correctFound = False

    for nrpsscore in candidateNRPSs_sorted:
        nrps = nrpsscore[0]
        score = nrpsscore[1]


