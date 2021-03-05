
# from read_gbk_files import *
import os
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re

def read_gbk_file_get_features(gb_file):
    from Bio import SeqIO
#     gb_file = "LOIC01000001.1.cluster017.gbk"
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        # now do something with the record
#         print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
#         print repr(gb_record.seq)
        return gb_record



def index_genbank_features(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
    return answer
def index_genbank_clusters_clusterID(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :

                clusterID = feature.qualifiers['note'][0]
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                answer[clusterID] = index
#                 for value in feature.qualifiers[qualifier] :
#                     if value in answer :
#                         print "WARNING - Duplicate key %s for %s features %i and %i" \
#                            % (value, feature_type, answer[value], index)
#                     else :
#                         answer[value] = index
    return answer
def index_nrps_genbank_features(gb_record, feature_type, qualifier) :
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    predictions = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                if "Type" in feature.qualifiers['sec_met']:
                        product_type = info.split(":")[1][1:]
                        if product_type !='nrps':
                            continue
                for value in feature.qualifiers['locus_tag'] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
                        predictions[value] = feature.qualifiers['aSProdPred'][0]
    return answer, predictions

#Functions to read NRPSpredictor2 codes

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


def readCodes(codesFileAddress,topAA,minScore,cyclominerPath): #we do not filter at this point on any point. JUST READ the codes file.
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashPredicts ={}
    antiSmashPredictsScore = {}
    antiSmashPredictsRank = {}
    import numpy as np
    gbkFile = open(codesFileAddress,"r")
    line =gbkFile.readline().strip()
#     new_contig = line.split()[0].split("_A")[0] #these are orfs not contigs
    new_contig = line.split()[0] #these are aDomains
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
        antiSmashPredictsRank[contig] = {}

        while new_contig ==contig:
            antiSmashPredicts[contig][pos] = []
            antiSmashPredictsScore[contig][pos] = []
            antiSmashPredictsRank[contig][pos] = []
            old_score =110
            rank=-1

            for i in range(15):
                # if  rank>topAA:
                #     print rank
                #     continue
                if len(set(antiSmashPredictsScore[contig][pos]))>topAA:
                    continue
                aminoScore = aminoAcids[i]
                aminoName = aminoScore.split("(")[0]
                score = float(aminoScore.split("(")[1].split(")")[0])
                if score<old_score:
                        rank+=1
                        if rank==3:
                            break
                if i==0:
                    adomainmaxscore = score
                if aminoName in aminFullName2Mass:
                    if score<minScore:
                        continue
                    antiSmashPredicts[contig][pos].append(aminFullName2Mass[aminoName])

                    newscore = int(round(float(score)/float(adomainmaxscore)*1.0,2)*100)
                    # if newscore < 80:
                    #     break
                    antiSmashPredictsScore[contig][pos].append(newscore)
                    # if score<old_score:
                    #     rank+=1
                    antiSmashPredictsRank[contig][pos].append(rank)
                old_score = score
            pos += 1
            line =gbkFile.readline().strip()
            if not line: 
                break
            lineSplit = line.split()
#             new_contig =  lineSplit[0].split("_A")[0]
            new_contig =  lineSplit[0]
            predictions = line.split()[2]
            aminoAcids = predictions.split(";")
    finalAntiSmashPredicts = {}
    finalAntiSmashPredictsScore = {}
    finalAntiSmashPredictsRank = {}
    for contig in antiSmashPredicts:
        if len(antiSmashPredicts[contig]) <0:
            continue
        finalAntiSmashPredicts[contig] = antiSmashPredicts[contig]
        finalAntiSmashPredictsScore[contig] = antiSmashPredictsScore[contig]
        finalAntiSmashPredictsRank[contig] = antiSmashPredictsRank[contig]
    return finalAntiSmashPredicts,finalAntiSmashPredictsScore,finalAntiSmashPredictsRank



def get_final_nrpPred_per_custer(gb_record,delete_orf):
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    #dictionary gets the  
    ctg_orf_Anum_aSDomain_index = index_genbank_features(gb_record,"aSDomain","label")
    # gets the CDs that has all the info regarding antismash-predicted peptide
    ctg_nrp_index, ctg_nrp_prediction = index_nrps_genbank_features(gb_record,"CDS","aSProdPred") 
    clusterName = gb_record.name
    allORF_NRPs_domain_names= {}
    allORF_NRPs_domain_locs= {}
    for cds in ctg_nrp_index: 
        allORF_NRPs_domain_names[cds] = []
        allORF_NRPs_domain_locs[cds] = []
        alldomains = []
        for info in gb_record.features[ctg_nrp_index[cds]].qualifiers['sec_met']:
            start = gb_record.features[ctg_nrp_index[cds]].location.nofuzzy_start
            end = gb_record.features[ctg_nrp_index[cds]].location.nofuzzy_end
            strand = gb_record.features[ctg_nrp_index[cds]].strand
            if "Type" in info:
                product_type = info.split(":")[1][1:]
            if 'NRPS/PKS Domain: ' in info:
                domainInfo = info.split()
                domain = domainInfo[2]
                if domain == 'AMP-binding':
                    domain= [a for a in domainInfo if 'ctg' in a][0][:-1]
                alldomains.append(domain)
        if gb_record.features[ctg_nrp_index[cds]].strand == -1:
            allORF_NRPs_domain_names[cds] = alldomains
            print alldomains
        else:
            allORF_NRPs_domain_names[cds] = alldomains
            # print alldomains
        allORF_NRPs_domain_locs[cds].append((start,end,strand))
    import re
    def atoi(text):
        return int(text) if text.isdigit() else text
    def natural_keys(text):
        return [ atoi(c) for c in re.split('(\d+)', text) ]

    # sorted_ctgs = allORF_NRPs_domain_locs
    sorted_ctgs = [x for x in allORF_NRPs_domain_locs.keys()]
    sorted_ctgs.sort(key=natural_keys) #sorted list of ORF names
    previous_loc = (-10000,-10000,-10000)
    finish_group = False
    new_group = False

    all_nrp_groups = {}
    all_nrp_groups_locs = {}
    group_num = -1
    for ctg in  sorted_ctgs:  # find the groups based on conditions new_group only if: 
                                #1) the location of previous domain is over 10000  bp away.
                                #2) the  TE domain was the previous domain
                                #3) the previous ORF was not on the same strand
        new_loc = allORF_NRPs_domain_locs[ctg][0]
        if finish_group == True:
            new_group= True
        else: 
            new_group = False
        finish_group = False
        if new_loc[2] != previous_loc[2]:
            new_group = True
        distance = new_loc[0]-previous_loc[1] 
        if not new_group and distance>1000:
            new_group = True
        if new_loc[2]  == 1:
            if allORF_NRPs_domain_names[ctg][-1] == 'Thioesterase':
#turn on for ending                finish_group = True
                finish_group = True
        elif allORF_NRPs_domain_names[ctg][0] == 'Thioesterase':
#turn on for ending                finish_group = True
            new_group = True

        previous_loc = new_loc

        if new_group:
            group_num += 1
            all_nrp_groups[group_num] = [ctg]
            all_nrp_groups_locs[group_num] =  new_loc
        else:
            all_nrp_groups[group_num].append(ctg)
            changing_last_loc = list(all_nrp_groups_locs[group_num])
            new_tuple = (all_nrp_groups_locs[group_num][0],new_loc[1],new_loc[2])
            all_nrp_groups_locs[group_num] = new_tuple
            if new_loc[2] != all_nrp_groups_locs[group_num][2]:
                print "There's something wrong with strand  information"

    final_nrps_domains = dict()
    import itertools
    totalnumCombs = 0 
    linear = False
    for grp in all_nrp_groups:
        #Mixing ORFs along the NRPS

        #choosing combinations of different sizes in the orfs
        if len(all_nrp_groups[grp]) >delete_orf:
            for comb in itertools.combinations(all_nrp_groups[grp], len(all_nrp_groups[grp])-delete_orf):
                totalnumCombs +=1 
                group = totalnumCombs
                final_nrps_domains[group] = []
                for ctg in comb:
                    final_nrps_domains[group] += allORF_NRPs_domain_names[ctg]

        # if len(all_nrp_groups[grp]) >2:
        #     for comb in itertools.combinations(all_nrp_groups[grp], len(all_nrp_groups[grp])-1):
        #         totalnumCombs +=1 
        #         group = totalnumCombs
        #         final_nrps_domains[group] = []
        #         for ctg in comb:
        #             final_nrps_domains[group] += allORF_NRPs_domain_names[ctg]
        #         if all_nrp_groups_locs[grp][2] == -1:
        #             final_nrps_domains[group].reverse()
        # totalnumCombs +=1 
        # group = totalnumCombs
        
        # for ctg in all_nrp_groups[grp]:
        #     final_nrps_domains[group] += allORF_NRPs_domain_names[ctg]
        # if all_nrp_groups_locs[grp][2] == -1:
        #     final_nrps_domains[group].reverse()
    return final_nrps_domains
        




        
def final_putative_NRP_gbk_code(
        adomain_predictions,
        adomain_predictions_scores,
        adomain_predictions_ranks,
        adomain_predictions_byrank,
        nrp_group,allSpecMasses,
        threshold,
        delete_aas,
        maxMod):
    import itertools
    final_nrps_info_local = {}
    max_mass_threshold = 1800
    for possibleMass in [round(p/100.0,2) for p in range(40000, max_mass_threshold*100)]:
        final_nrps_info_local[round(possibleMass,2)] = {}
    from collections import deque

    def findCyclicPermutation(a):
        return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
    
    final_num_nrps = 0

    def checkMass(nrpseq,nrp_name):
        #checks if the putative nrp matches the mass condition and a spectrum 
        mass = round(sum(nrpseq),2)
        if mass >400 and mass <( max_mass_threshold):
            if allSpecMasses[mass] == 1:
                final_nrps_info_local[mass].update({nrp_name:tuple(findCyclicPermutation((nrpseq)))})
    import itertools

    def create_nrp_seq_name(structure, ctgNum, groupName, compound_number):
        #creates group name and the NRP name: the NRP name is a mixture 
        nrp_name = ctgNum + "_" + groupName + "_" + str(compound_number).zfill(6)
        nrp_tuple = tuple([j for j in (n for n ,_ in structure)])
        return nrp_name, nrp_tuple

    def checkSizeofCoreDB(adomain_predictions_grp,grp,i): 
        #counts the maximum number of putative NRPs 
        product = 1
                # print groupName
        for x in [len([a for (a,b,c) in adomain_predictions_grp[orf] if c<=i]) for orf in grp]:
            product *= x
        return product
    def _generate_combination_of_rank(maximumrank, listofranks, length):
        def crossprod (list1, list2):
            output = 0
            for i in range(0,len(list1)):
                output += list1[i]*list2[i]

            return output

        def breakit(target, coins):
            final_targets = []
            coinslimit = [(target / coins[i]) for i in range(0,len(coins))]
            count = 0
            temp = []
            for i in range(0,len(coins)):
                temp.append([j for j in range(0,coinslimit[i]+1)])
            r=[[]]
            for x in temp:
                t = []
                for y in x:
                    for i in r:
                        t.append(i+[y])
                r = t

            for targets in r:
                if crossprod(targets, coins) == target:

                    final_targets.append(targets)
                    count +=1

            return count, final_targets
        final_rank_combinations = {}
        final_num_comb = 0
        import itertools
        allpermutations = set()
        coins = [1,2]
        totalrank = maximumrank
        for rank in range(totalrank + 1): #this shoudl bethe length of the structure
            num_comb, final_rank_combinations[rank] = breakit(rank, coins)
            final_num_comb += num_comb
            for combination in final_rank_combinations[rank]:
                num_zeros = length - sum(combination)
                combination_listof_ranks= [0]*num_zeros + [1]*combination[0] + [2]*combination[1]

                for one_positions in itertools.combinations(range(length),combination[0]):
                    leftpositions = range(length)[:]
                    [leftpositions.remove(p) for p in list(one_positions)]

                    for two_positions in itertools.combinations(leftpositions,combination[1]):
                        permutation = [0]*length
                        for pos in one_positions:
                            permutation[pos] = 1 
                        for pos in two_positions:
                            permutation[pos] = 2
                    allpermutations.add(tuple(permutation))
        return allpermutations

    def delet_Adomains_create_nrps(NRPS_adomains, delete_aas, ctgNum):        
        #deletes delete_amino number of A-domains from the NRPS
        compound_number = 0
        import numpy
        for i in range(delete_aas,delete_aas+1):
            if len(NRPS_adomains)<=i:
                continue

            for group in itertools.combinations(NRPS_adomains, len(NRPS_adomains)-i):
                #name of the group is the contig name succeeded by the name of orfs and A-domains involved
                groupName = "-".join(tuple(["orf"+adomain.split("orf")[1] for adomain in group]))
                compound_number_grp = 0
                nrp_rank_number ={}
                product0 = 1
                max_bruteforce_count = {}
                maxrankfound = False
                for i in range(3):
                    max_bruteforce_count[i]= checkSizeofCoreDB(adomain_predictions_ranks, group,i)
                    # fins th4e smallest rank that can be considered where the maximum number of 
                    if max_bruteforce_count[i] > 10000 and not maxrankfound:
                        max_aa_rank = i
                        maxrankfound = True
                # largest_bruteforce_possible = max(max_bruteforce_count.values()) 
                if max_bruteforce_count[2]>1000000:
                    maxTotalRank = 3
                elif max_bruteforce_count[1]>10000:
                    maxTotalRank = 5 
                else:
                    maxTotalRank = 10
                if maxMod == 0:
                    maxTotalRank = len(NRPS_adomains)
                if not maxrankfound:
                    max_aa_rank = 0
                print "BRUTE FORCE NUMBER\t1\t" + str(max_bruteforce_count)
                scores = {}

                print "Finding the core NRPs ..."
                for vv in range(0,1):
                    finalnum = 0
                    for rank_permutation in _generate_combination_of_rank(maxTotalRank,range(min(max_aa_rank,1),max_aa_rank+1),len(group)):
                        # print "------"
                        # print rank_permutation 
                        feasible = True
                        totalnum = 0
                        product = 1
                        for i in range(len(rank_permutation)):
                            rank = rank_permutation[i]
                            if group[i] not in adomain_predictions_byrank[rank]:
                                feasible = False
                                break
                        if not feasible:
                            continue

                        for x in [len([a for (a,b,c) in adomain_predictions_byrank[rank_permutation[i]][group[i]] ]) for i in range(len(group))] :
                            product *= x
                        finalnum += product
                        totalnum = 0
                        for x in itertools.product(*[ [(a,b) for (a,b,c) in adomain_predictions_byrank[rank_permutation[i]][group[i]] ] for i in range(len(group)) ] ):
                            structure = tuple(x)
                            totalnum +=1 


                            #check if the score of amino acids meet the criteria
                            rankssum = round(round(sum(z for _, z in structure))/float(len(structure)),2)

                            if len(structure)>=11:
                                scorethresh = 95
                            elif len(structure)<11 and len(structure)>6:
                                scorethresh = 93
                            else:
                                scorethresh = 90
                            threshold=scorethresh

                            if rankssum >= threshold: 
                                nrp_name, nrp_tuple = create_nrp_seq_name(structure, ctgNum, groupName+"_"+str(rankssum), compound_number)
                                # print nrp_name
                                nrp_name += str(rankssum)
                                
                                nrptuple_str = "_" + "-".join([str(str(int(a))) for a in nrp_tuple])
                                nrp_name += nrptuple_str
                                checkMass(nrp_tuple, nrp_name)
                                compound_number +=1
                                compound_number_grp +=1
                    # print "BruteForce-vs-NRPminer\t" + groupName + "\t" + str(product) + "\t" + str(compound_number_grp)
                    # threshold = sorted(scores.values(), reverse = True)[min(1000,len(scores)-1)]
                    # threshold=scorethresh
                    # if rankssum >= threshold: 
                    # for structure in [key for (key,value) in sorted(scores.items(), reverse = True)[0:min(1000,len(scores)-1)]]:
                    #     rankssum = scores[structure]
                    
                    #     nrp_name, nrp_tuple = create_nrp_seq_name(structure, ctgNum, groupName+"_"+str(rankssum), compound_number)
                    #     # print nrp_name
                    #     nrp_name += str(rankssum)
                    #     checkMass(nrp_tuple, nrp_name)
                    #     compound_number +=1
                    #     compound_number_grp +=1
                    print "-----"
                    print group
                    print "BruteForce-vs-NRPminer\t" + groupName +"\t" + "\t".join([str(max_bruteforce_count[z]) for z in max_bruteforce_count]) +"\t" + str(compound_number_grp) + "\t" + str(maxTotalRank) + "\t" + str(max_aa_rank)
                # print sorted(scores.items(), key=lambda x: x[1])
        # exit()
        return compound_number
    def add_Adomains_create_nrps(NRPS_adomains, delete_aas, ctgNum):
        #deletes delete_amino number of A-domains from the NRPS
        compound_number = 0
        print "==================================="
        print "==================================="
        # print NRPS_adomains
        for i in range(delete_aas,delete_aas+1):
            if len(NRPS_adomains)<=i:
                continue
            for group in itertools.combinations(NRPS_adomains, len(NRPS_adomains)-i):
                #name of the group is the contig name succeeded by the name of orfs and A-domains involved
                groupName = "-".join(tuple(["orf"+adomain.split("orf")[1] for adomain in group]))
                # print [len(adomain_predictions[orf]) for orf in group]
                product = 1
                for x in [len(adomain_predictions[orf]) for orf in group]:
                    product *= x
                # print product
                for x in itertools.product(*[adomain_predictions[orf] for orf in group ] ):
                    structure = tuple(x)
                    # print structure
                    #check if the number of amino acids with rank 1 and 2 doesn't get beyond the thershold
                    rankssum = sum(z for _, z in structure)/float(len(structure))                    
                    # rankssum = sum(z for _, z in structure)
                    if rankssum <= threshold:
                        nrp_name, nrp_tuple = create_nrp_seq_name(structure, ctgNum, groupName+"_"+str(rankssum), compound_number)
                        nrp_name += str(rankssum)
                        checkMass(nrp_tuple, nrp_name)
                        compound_number +=1
        return compound_number
    # maxRank = 10
    for NRPS_adomains in nrp_group:
        ctgNum = NRPS_adomains[0].split("_")[0]
        compound_counts = delet_Adomains_create_nrps(NRPS_adomains,delete_aas,ctgNum)
    return final_nrps_info_local


def create_putative_nrp_dataset(antismash_res_dir, allSpecMasses,cyclominerPath,threshold,delete_orf,delete_aas,topaa,minscore,maxMod):
    # antismash_res_dir = "/Users/bahar/workspace/npd_tools/cyclominer_bloom/antismash_results/xeno_antismash/LOIC01.1.fasta_out"
    # antismash_res_dir = "/Users/bahar/workspace/npd_tools/cyclominer_bloom/antismash_results/xeno_antismash/NIUA01.fasta_out"
    import glob
    all_gbk_files_in_dir = glob.glob(antismash_res_dir+"/*cluster*.gbk")

    per_ctgs_NRP_groups_domains = dict()
    per_ctgs_code_File = dict()
    def get_final_nrp_groups_per_ctg(per_cluster_final_nrps_domains):
        per_ctgs_code_File = dict()
        total_num_groups = 0
        for group in per_cluster_final_nrps_domains:
            domains = per_cluster_final_nrps_domains[group]
            PKS_domains = [domain for domain in domains if re.match(r'PKS.*', domain)]
            if len(PKS_domains)>3:
                continue
            all_Adomains = [domain for domain in domains if re.match(r'ctg.*orf.*', domain)]
            numAA = len(all_Adomains)
            if numAA<3 or numAA>16:
                continue
            total_num_groups += 1
            ctg_num = all_Adomains[0].split("_")[0]
            #if ctg_num != 'ctg35':
            #    continue
            if ctg_num not in per_ctgs_NRP_groups_domains:
                per_ctgs_NRP_groups_domains[ctg_num] = [tuple(all_Adomains)]
            else:
                per_ctgs_NRP_groups_domains[ctg_num].append(tuple(all_Adomains))

        return per_ctgs_NRP_groups_domains, per_ctgs_code_File, total_num_groups

    for gb_file in all_gbk_files_in_dir:
        print gb_file
        gb_record = read_gbk_file_get_features(gb_file)
        clusterName = gb_record.name
        per_cluster_final_nrps_domains = get_final_nrpPred_per_custer(gb_record,delete_orf)
        test, per_ctgs_code_File, total_num_groups = get_final_nrp_groups_per_ctg(per_cluster_final_nrps_domains)

    all_nrpspredictorCode_files_in_dir = glob.glob(antismash_res_dir+"/nrpspks_predictions_txt/ctg*_nrpspredictor2_codes.txt")
    per_ctg_nrpspred_codes = {}
    per_ctg_nrpspred_scores = {}
    per_ctg_nrpspred_ranks = {}
    per_ctg_nrpspred_byrank = {}
    for code_file in all_nrpspredictorCode_files_in_dir:
        per_orf_nrpspredictor2_codes,per_orf_nrpspredictor2_scores,per_orf_nrpspredictor2_ranks =  readCodes(code_file,topaa,minscore,cyclominerPath)

        for orf in per_orf_nrpspredictor2_codes:
            ctg = orf.split("_")[0]

            if ctg in per_ctg_nrpspred_scores:
                # per_ctg_nrpspred_codes[ctg].update({orf: [(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i]) for i in range(len(per_orf_nrpspredictor2_codes[orf][0]))] })
                per_ctg_nrpspred_scores[ctg].update({orf:per_orf_nrpspredictor2_scores[orf]})

                per_ctg_nrpspred_ranks[ctg].update({orf: [(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i]) for i in range(len(per_orf_nrpspredictor2_codes[orf][0]))] })
                for i in range(len(per_orf_nrpspredictor2_codes[orf][0])):
                    current_rank = per_orf_nrpspredictor2_ranks[orf][0][i]

                    if orf in per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]]:
                        per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]][orf].append((per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i]))
                    else:
                        per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]].update({orf:[(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i])]})
                    # per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]].update({orf:[(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i]) for i in range(len(per_orf_nrpspredictor2_codes[orf][0]))]})

            else:
                per_ctg_nrpspred_codes[ctg] = {orf: [(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i]) for i in range(len(per_orf_nrpspredictor2_codes[orf][0]))] }
                per_ctg_nrpspred_scores[ctg] = {orf:per_orf_nrpspredictor2_scores[orf]}
                per_ctg_nrpspred_ranks[ctg] = {orf:[(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i]) for i in range(len(per_orf_nrpspredictor2_codes[orf][0]))]}
                if ctg not in per_ctg_nrpspred_byrank:  
                    per_ctg_nrpspred_byrank[ctg] = {}
                    # per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]] = {}
                for i in range(3):
                    per_ctg_nrpspred_byrank[ctg][i] = {}
                for i in range(len(per_orf_nrpspredictor2_codes[orf][0])):
                    current_rank = per_orf_nrpspredictor2_ranks[orf][0][i]
                    per_ctg_nrpspred_byrank[ctg][per_orf_nrpspredictor2_ranks[orf][0][i]].update({orf:[(per_orf_nrpspredictor2_codes[orf][0][i],per_orf_nrpspredictor2_scores[orf][0][i],per_orf_nrpspredictor2_ranks[orf][0][i])]})


    final_nrps_info = dict()
    max_mass_threshold = 1800
    for possibleMass in [round(p/100.0,2) for p in range(40000, max_mass_threshold*100)]:
        final_nrps_info[round(possibleMass,2)] = {}
    for ctg in per_ctg_nrpspred_codes:
        if ctg in per_ctgs_NRP_groups_domains:
            ctg_nrps_info = final_putative_NRP_gbk_code(
                per_ctg_nrpspred_codes[ctg],per_ctg_nrpspred_scores[ctg],per_ctg_nrpspred_ranks[ctg],per_ctg_nrpspred_byrank[ctg],
                per_ctgs_NRP_groups_domains[ctg],allSpecMasses,threshold,delete_aas,maxMod)
            [final_nrps_info[mass].update(ctg_nrps_info[mass]) for mass in ctg_nrps_info]
    return final_nrps_info
# create_putative_nrp_dataset("/Users/bahar/workspace/npd_tools/cyclominer_bloom/antismash_results/xeno_antismash/LOIC01.1.fasta_out", "")

def output_nrpsgraph(candidateNRPSs, output):
    #This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
    from operator import itemgetter
    def createPNPGraphLines(pnp,name):

        num = len(pnp)
        lines= []
        nameLine = name
        lines.append(nameLine) 
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

    # def createPNPGraphLinesLinear(pnp,name):

    #     num = len(pnp)
    #     lines= []
    #     nameLine = name+"linear"
    #     lines.append(nameLine)
    #     newline = "number of components : "+ str(num)
    #     lines.append(newline)
    #     for i in range(num):
    #         aa = pnp[i]
    #         newline = str(i) + " CXHXNX " + str(aa)
    #         lines.append(newline)
    #     newline = "number of bonds : " + str(num-1)
    #     lines.append(newline)
    #     for i in range(num-1):
    #         newline = str(i) + " -NC> " + str(i+1)
    #         lines.append(newline)
    #     #newline = str(num-1) + " -NC> " + str(0)
    #     lines.append(newline)
    #     return lines

    # candidateNRPSs_sorted = sorted(candidateNRPSs.items(), key= itemgetter(1), reverse=True)

    scores = set()
    
    peptidesUptoCorrectScore = []
    correctFound = False
    foundOne = False
    for nrpPair in candidateNRPSs:
        nrps = nrpPair[1]
        compoundName = nrpPair[0]
        # score = nrpsscore[1]
        # score = 3
        # if score>2:
        graphFile = open(output+"_candidateNRPs.GRAPHS.txt","a")
        #'''
        linesToAdd = createPNPGraphLines(nrps,compoundName)
        foundOne = True
        for line in linesToAdd:
            graphFile.write(line+"\n")
        #'''linesToAdd = []
        # linesToAdd = createPNPGraphLinesLinear(nrps,compoundName)
        for line in linesToAdd:
            graphFile.write(line+"\n")

    return foundOne
