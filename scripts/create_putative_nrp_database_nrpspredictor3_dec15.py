
# from read_gbk_files import *
import os
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re

#function to process the gbk files
def read_gbk_file_get_features(gb_file):
    '''read the GBK files'''
    from Bio import SeqIO
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        return gb_record


def index_genbank_features(gb_record, feature_type, qualifier) :
    #function to retrieve the gbk features 
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



def get_specificty_genbank_features(
    gb_record,
    predictor='Stachelhaus') :
    # fiven a gbk record finds the specificity prediction based on predictor 
    # predictor between ['Stachelhaus','NRPSpredictor3','pHMM','SANDPUMA'])
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = ['no_call']
    feature_type = 'aSDomain'
    domain = 'AMP-binding'
    predictions = {}

    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type : #type: aSDomain
            if "domain" in feature.qualifiers:
                if feature.qualifiers['domain'][0] == domain:
                    if "specificity" in feature.qualifiers:
                        for prediction in feature.qualifiers["specificity"]:
                            if predictor in prediction:
                                if 'label' in feature.qualifiers:
                                    Adomain_ID= feature.qualifiers['label'][0]
                                    answer = prediction.split(":")[1][1:]
                                    predictions[Adomain_ID] = answer.split("|")
    return predictions


def index_genbank_clusters_clusterID(
    gb_record,
    feature_type,
    qualifier
    ):
    #get the features based on the BGC cluster
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
                answer[clusterID] = index

    return answer

def index_nrps_genbank_features(
    gb_record,
    feature_type,
    qualifier) : #finds different assembly lines
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = dict()
    predictions = dict()
    print qualifier

    feature_type = 'aSDomain'
    for (index, feature) in enumerate(gb_record.features) : #go through the features in gbk 
        if feature.type=="aSDomain":
            print index
            print feature
            print feature.type
            # print feature.qualifiers
            print feature.qualifiers['aSDomain']
            
            # for x in feature.qualifiers:
            #     print "--"
            #     print x
            #     print feature.qualifiers[x][0]
            # exit()
        # print feature.qualifiers['aSDomain']
        # exit()
        if feature.type==feature_type :
            print feature.qualifiers

            if 'aSDomain' in feature.qualifiers:
                print "HELLLO?"
            # if qualifier in feature.qualifiers :
                print feature.qualifiers

                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries

                domain_id = feature.qualifiers['domain_id'][0]
                # genome_start = feature.qualifiers['domain_id']
                # end = feature.qualifiers['domain_id']
                # strand = 

                location = feature.location
                
                start= feature.location.nofuzzy_start
                end = feature.location.nofuzzy_end
                strand = feature.strand

                
                # if "Type" in feature.qualifiers['sec_met']:
                #         product_type = info.split(":")[1][1:]
                #         if product_type !='nrps':
                #             continue
                for value in feature.qualifiers['locus_tag'] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        #answer includes all NRPS domains with start, end, strand and then all qualifier info
                        answer[domain_id] = (start, end, strand, feature.qualifiers, feature.qualifiers['locus_tag'][0])
                        # answer[value] = index
                        # predictions[value] = feature.qualifiers['aSProdPred'][0]
                        predictions = []

    return answer, predictions





def get_final_nrpPred_per_custer(gb_record,delete_orf): 
    # readGBK
    #read GBK record and find NRPS assembly lines
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re

    ctg_orf_Anum_aSDomain_index = index_genbank_features(gb_record,"aSDomain","domain_id")
    # print ctg_orf_Anum_aSDomain_index

    ctg_nrpdomains_info , ctg_nrp_prediction = index_nrps_genbank_features(gb_record,"aSDomain","domain_id")
    print "+++++++++"
    print ctg_nrpdomains_info
    # exit()
    clusterName = gb_record.name
    allORF_NRPs_domain_names= {}
    allORF_NRPs_domain_locs= {}
    allORF_NRPs_domain_predictions = {}
    def Reverse(lst): 
        return [ele for ele in reversed(lst)] 
    strand = 1
    for domain in ctg_nrpdomains_info:  #cds is the domain_id
        
        start,end,strand,info,cds = ctg_nrpdomains_info[domain]
        # alldomains.append(domain)
        if cds not in allORF_NRPs_domain_names:
            allORF_NRPs_domain_names[cds] = []
            allORF_NRPs_domain_locs[cds] = []
            allORF_NRPs_domain_predictions[cds] = [] #tuples with (stach,SVM,pHMM,SANDPUM)
            alldomains = []
        
        allORF_NRPs_domain_names[cds].append(domain)
        allORF_NRPs_domain_locs[cds].append((start,end,strand))

        #         if domain == 'AMP-binding':
        #             domain= [a for a in domainInfo if 'ctg' in a][0][:-1]
        #             # if "Substrate specificity predictions" in info:
        #             #     stach_code = info
        #             #     domain= [a for a in domainInfo if 'ctg' in a][0][:-1]
        #         alldomains.append(domain)


        # print gb_record.features[ctg_nrp_index[cds]]
        # # exit()


        #     if 'NRPS/PKS Domain: ' in info:
        #         domainInfo = info.split()
        #         domain = domainInfo[2]
        #         if domain == 'AMP-binding':
        #             domain= [a for a in domainInfo if 'ctg' in a][0][:-1]
        #             # if "Substrate specificity predictions" in info:
        #             #     stach_code = info
        #             #     domain= [a for a in domainInfo if 'ctg' in a][0][:-1]
        #         alldomains.append(domain)
        
        # '''
    for cds in allORF_NRPs_domain_locs:
        print allORF_NRPs_domain_locs
        if allORF_NRPs_domain_locs[cds][-1] == -1:
            listofdomains= allORF_NRPs_domain_names[cds]
            allORF_NRPs_domain_names[cds] = listofdomains[::-1]
        # else:
        #     allORF_NRPs_domain_names[cds] = alldomains
    print allORF_NRPs_domain_names[cds]

    
    
    exit()
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
                finish_group = True
        elif allORF_NRPs_domain_names[ctg][0] == 'Thioesterase':
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


    return final_nrps_domains
        



def final_putative_NRP_gbk_code(
        adomain_predictions_byrank,
        nrp_group,allSpecMasses,
        threshold,
        delete_aas,
        maxMod):
    import itertools
    final_nrps_info_local = {}
    max_mass_threshold = 1800
    for possibleMass in [round(p/100.0,2) for p in range(30000, max_mass_threshold*100)]:
        final_nrps_info_local[round(possibleMass,2)] = {}
    from collections import deque

    def findCyclicPermutation(a):
        return sorted([ a[n:] + a[:n] for n in range(len(a)) ])[0]
    
    final_num_nrps = 0

    def checkMass(nrpseq,nrp_name):
        #checks if the putative nrp matches the mass condition and a spectrum 
        mass = round(sum(nrpseq),2)

        if mass > 300 and mass <( max_mass_threshold):
            if allSpecMasses[mass] == 1:
                final_nrps_info_local[mass].update({nrp_name:tuple(findCyclicPermutation((nrpseq)))})
    import itertools

    def create_nrp_seq_name(structure, ctgNum, groupName, compound_number):
        #creates group name and the NRP name: the NRP name is a mixture 
        nrp_name = ctgNum + "_" + groupName + "_" + str(compound_number).zfill(6)
        nrp_tuple_name = tuple([j for j in (n for _ ,_ , n in structure)])
        nrp_tuple_mass = tuple([j for j in (n for n ,_ , _ in structure)])        
        return nrp_name, nrp_tuple_name, nrp_tuple_mass

    def checkSizeofCoreDB(adomain_predictions_grp,grp,i): 
        #counts the maximum number of putative NRPs 
        product = 1
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
                groupName = "-".join(tuple([adomain.split("_")[1]+"_"+adomain.split("_")[2] for adomain in group]))
                product = 1
                for i in range(len(group)):
                    total_num_grp = 0
                    for r in adomain_predictions_byrank:
                        if group[i] in adomain_predictions_byrank[r]:
                            total_num_grp += len([a for (a,b,c,d) in adomain_predictions_byrank[r][group[i]] ])
                    product *= total_num_grp                
                #print "BRUTE FORCE NUBMER is " + str(product) + "\t" + groupName
                bruteforce_num = product
 
                compound_number_grp = 0
                nrp_rank_number ={}
                product0 = 1
                max_bruteforce_count = {}
                maxrankfound = False
                if not maxrankfound:
                    max_aa_rank = 0
                
                maxTotalRank = 0
                scores = {}
                finalnum = 0

                for totalRank in range(len(group)): 
                    if finalnum > 10000:
                        break
                    for rank_permutation in _generate_combination_of_rank(totalRank,range(min(max_aa_rank,1),max_aa_rank+1),len(group)):
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
                        product = 1
                        for x in [len([a for (a,b,c,d) in adomain_predictions_byrank[rank_permutation[i]][group[i]] ]) for i in range(len(group))] :
                            product *= x
                        finalnum += product
                        if finalnum > 100000:
                            if finalnum> 10000:
                                maxTotalRank = totalRank-1
                            else:
                                maxTotalRank = totalRank
                            break
                        maxTotalRank = totalRank

                print "Finding the core NRPs ..."
                # maxTotalRank = 20
                for vv in range(0,1):
                    finalnum = 0
                    add2onehundred = False
                    min100score = 0 
                    firstonehundred = {}
                    firstonehundred[min100score] = []
                    for rank_permutation in _generate_combination_of_rank(maxTotalRank,range(min(max_aa_rank,1),max_aa_rank+1),len(group)):
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

                        for x in [len([a for (a,b,c,d) in adomain_predictions_byrank[rank_permutation[i]][group[i]] ]) for i in range(len(group))] :
                            product *= x
                        finalnum += product
                        totalnum = 0
                        for x in itertools.product(*[ [(a,b,d) for (a,b,c,d) in adomain_predictions_byrank[rank_permutation[i]][group[i]] ] for i in range(len(group)) ] ):
                            structure = tuple(x)
                            totalnum +=1 
                            add2onehundred = False
                            # rankssum = round(round(sum(z for _, z in structure))/float(len(structure)),2)
                            rankssum = round(round(sum(z for _, z, _ in structure)),0)
                            if min100score < rankssum:
                                add2onehundred = True
                                if len(firstonehundred) >= 1000:
                                    del firstonehundred[min100score]

                            if add2onehundred:
                                if rankssum not in firstonehundred:
                                    firstonehundred[rankssum] = []

                                firstonehundred[rankssum].append(structure)
                                min100score = min(firstonehundred.keys())

                    total = 0
                    # minscore_to_check = 90
                    maxscoreobserved= max(firstonehundred.keys())
                    minscore_to_check = 0.6*maxscoreobserved
                    for k in sorted(firstonehundred.keys(),reverse=True):
                        total += len(firstonehundred[k])
                        if total >1000:
                            minscore_to_check  = k
                            break
                    for score in sorted(firstonehundred,reverse=True):
                            if score < minscore_to_check:
                                break
                            for structure in firstonehundred[score]:
                                nrp_name, nrp_tuple_name, nrp_tuple_mass = create_nrp_seq_name(structure, ctgNum, groupName+"_"+str(score), compound_number)
                                nrp_list = list(nrp_tuple_mass)
                                nrptuple_name_str = "_" + "-".join([str(str((a))) for a in nrp_tuple_name])
                                nrp_name += nrptuple_name_str
                                checkMass(nrp_tuple_mass, nrp_name)
                                compound_number +=1
                                compound_number_grp +=1
                    print "NRPMINER NUMBER\t:" + str(group) + "\t" + str(total) + "\t" + str(minscore_to_check)
        return compound_number

    for NRPS_adomains in nrp_group:
        ctgNum = NRPS_adomains[0].split("_")[0]
        compound_counts = delet_Adomains_create_nrps(NRPS_adomains,delete_aas,ctgNum)
    
    return final_nrps_info_local


def create_putative_nrp_dataset(antismash_res_dir
    ,allSpecMasses
    ,cyclominerPath
    ,threshold
    ,delete_orf
    ,delete_aas
    ,topaa
    ,minscore
    ,maxMod):

    import glob
    all_gbk_files_in_dir = glob.glob(antismash_res_dir+"/*34*region*.gbk")
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
            all_Adomains = [domain for domain in domains if re.match(r'ctg.*A.*', domain)]
            numAA = len(all_Adomains)
            #minimum length of assembly lines considered is 3 and maximum is 20
            if numAA<3 or numAA>20:
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
    ## ERRROR HERE
    print all_gbk_files_in_dir
    for gb_file in all_gbk_files_in_dir:
        gb_record = read_gbk_file_get_features(gb_file)
        clusterName = gb_record.name
        per_cluster_final_nrps_domains = get_final_nrpPred_per_custer(gb_record,delete_orf)
        print per_cluster_final_nrps_domains
        test, per_ctgs_code_File, total_num_groups = get_final_nrp_groups_per_ctg(
            per_cluster_final_nrps_domains)


    all_nrpspredictorCode_files_in_dir = glob.glob(antismash_res_dir+"/nrpspks_predictions_txt/ctg*_nrpspredictor2_codes.txt")
    all_nrpspredictorSVM_files_in_dir = glob.glob(antismash_res_dir+"/nrpspks_predictions_txt/ctg*_nrpspredictor*_svm.txt")    
    all_nrpspredictorInd_files_in_dir = glob.glob(antismash_res_dir+"/nrpspks_predictions_txt/ctg*_ind.res.tsv")
    
    
    per_ctg_nrpspred_codes = {}
    per_ctg_nrpspred_scores = {}
    per_ctg_nrpspred_ranks = {}
    per_ctg_nrpspred_byrank = {}
    per_orf_nrpspredictor2_codes ={}
    per_orf_nrpspredictor2_ ={}
    per_orf_nrpspredictor2_scores = {}
    per_orf_nrpspredictor2_ranks = {}
    per_orf_nrpspredictor2_AAname = {}
    from database_functions import readCodes
    from database_functions import readSVMs
    from database_functions import readINDs
    from database_functions import readGBK_predictions

    for code_file in all_nrpspredictorCode_files_in_dir:
        per_orf_nrpspredictor2_codes_toadd,per_orf_nrpspredictor2_scores_toadd,per_orf_nrpspredictor2_ranks_toadd,per_orf_nrpspredictor2_AAname_toadd =  readCodes(code_file,topaa,minscore,cyclominerPath)
        if len(per_orf_nrpspredictor2_codes_toadd):            
            per_orf_nrpspredictor2_codes.update(per_orf_nrpspredictor2_codes_toadd)
            per_orf_nrpspredictor2_scores.update(per_orf_nrpspredictor2_scores_toadd)
            per_orf_nrpspredictor2_ranks.update(per_orf_nrpspredictor2_ranks_toadd)
            per_orf_nrpspredictor2_AAname.update(per_orf_nrpspredictor2_AAname_toadd)


    # normalizign the specificty scores.
    global per_orf_nrpspredictor2_SVMs
    global per_orf_nrpspredictor2_scores_SVM
    global per_orf_nrpspredictor2_AAname_SVM
    per_orf_nrpspredictor2_SVMs = {}
    per_orf_nrpspredictor2_scores_SVM = {}    
    per_orf_nrpspredictor2_AAname_SVM = {}
    ##################################
    ###### read SVM predictions from nrpspredictor*svm files under 
    ##################################
    for svm_file in all_nrpspredictorSVM_files_in_dir:
        per_orf_nrpspredictor2_SVMs_toadd,per_orf_nrpspredictor2_scores_SVM_toadd,per_orf_nrpspredictor2_ranks_SVM_toadd, per_orf_nrpspredictor2_AAname_SVM_toadd =  readSVMs(svm_file,topaa,minscore,cyclominerPath,antismash_res_dir)
        if len(per_orf_nrpspredictor2_SVMs_toadd): 
            per_orf_nrpspredictor2_SVMs.update(per_orf_nrpspredictor2_SVMs_toadd)
            per_orf_nrpspredictor2_scores_SVM.update(per_orf_nrpspredictor2_scores_SVM_toadd)
            per_orf_nrpspredictor2_AAname_SVM.update(per_orf_nrpspredictor2_AAname_SVM_toadd)


    ##### add IND file information in nrpspks
    for ind_file in all_nrpspredictorInd_files_in_dir: #this is to add to the NRPSpredictro3 SVM files. We add the other predictions from SANDPUMA and pHMM results.
        per_orf_nrpspredictor2_INDs_toadd,per_orf_nrpspredictor2_scores_INDs_toadd,per_orf_nrpspredictor2_ranks_INDs_toadd, per_orf_nrpspredictor2_AAname_INDs_toadd =  readINDs(ind_file,topaa,minscore,cyclominerPath)

        if len(per_orf_nrpspredictor2_INDs_toadd):
            for Adomain in per_orf_nrpspredictor2_INDs_toadd:
                if Adomain in per_orf_nrpspredictor2_SVMs:
                    for i in range(len(per_orf_nrpspredictor2_INDs_toadd[Adomain])):
                        amino = per_orf_nrpspredictor2_INDs_toadd[Adomain][i]
                        aminoName = per_orf_nrpspredictor2_AAname_INDs_toadd[Adomain][i]
                        if amino not in per_orf_nrpspredictor2_SVMs[Adomain]:
                            per_orf_nrpspredictor2_SVMs[Adomain].append(amino)
                            per_orf_nrpspredictor2_scores_SVM[Adomain].append(100)
                            per_orf_nrpspredictor2_AAname_SVM[Adomain].append(aminoName)
                        else:
                            ind = per_orf_nrpspredictor2_SVMs[Adomain].index(amino)
                            per_orf_nrpspredictor2_scores_SVM[Adomain][ind] = 100
                            per_orf_nrpspredictor2_AAname_SVM[Adomain][ind] = aminoName
                else:
                    per_orf_nrpspredictor2_SVMs[Adomain] = per_orf_nrpspredictor2_INDs_toadd[Adomain]
                    per_orf_nrpspredictor2_scores_SVM[Adomain] = per_orf_nrpspredictor2_scores_INDs_toadd[Adomain]
                    per_orf_nrpspredictor2_AAname_SVM[Adomain] = per_orf_nrpspredictor2_AAname_INDs_toadd[Adomain]

    ##################################
    ###### add GBK-read predictions.
    ##################################
    from database_functions import getBuildingBlocksForNRPSprediction    
    global aminFullName2Mass
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    from database_functions import nrpspredictor_clusters
    small_clusters_dict, large_clusters_dict = nrpspredictor_clusters()
    
    all_gbk_files_in_dir = glob.glob(antismash_res_dir+"/*34*region*.gbk")
    print all_gbk_files_in_dir
    def addAmino2dict(Adomain,aminoName,score):
        if aminoName == 'no_call':
            return
        if aminoName in aminFullName2Mass:
            amino = aminFullName2Mass[aminoName]
        else:
            return

        if Adomain not in per_orf_nrpspredictor2_SVMs:
            per_orf_nrpspredictor2_SVMs[Adomain] = []
            per_orf_nrpspredictor2_scores_SVM[Adomain] = []
            per_orf_nrpspredictor2_AAname_SVM[Adomain] = []

        if amino not in per_orf_nrpspredictor2_SVMs[Adomain]:
            per_orf_nrpspredictor2_SVMs[Adomain].append(amino)
            per_orf_nrpspredictor2_scores_SVM[Adomain].append(score)
            per_orf_nrpspredictor2_AAname_SVM[Adomain].append(aminoName)
        else:
            ind = per_orf_nrpspredictor2_SVMs[Adomain].index(amino)
            if per_orf_nrpspredictor2_scores_SVM[Adomain][ind] < score:
                per_orf_nrpspredictor2_scores_SVM[Adomain][ind] = score
            per_orf_nrpspredictor2_AAname_SVM[Adomain][ind] = aminoName
    
    for gbk_file in all_gbk_files_in_dir: #iterate through gbk files
        gb_record = read_gbk_file_get_features(gbk_file)
        for predictor in ['Stachelhaus','NRPSpredictor3','pHMM','SANDPUMA']:
            per_orf_gbk_predictions =  get_specificty_genbank_features(gb_record,'Stachelhaus')
            observed_small_clusters = []
            observed_large_clusters = []
            if len(per_orf_gbk_predictions):
                for Adomain in per_orf_gbk_predictions:
                    for aminoName in per_orf_gbk_predictions[Adomain]:
                        if aminoName in large_clusters_dict:
                            for lc_amino in large_clusters_dict[aminoName]: 
                                #lc_amino are the amino acids sharing the same cluster with aminoName
                                score = 80
                                addAmino2dict(Adomain,lc_amino,score)
                        if aminoName in small_clusters_dict:
                            for sc_amino in small_clusters_dict[aminoName]: 
                                score = 90
                                addAmino2dict(Adomain,sc_amino,score)
                        score = 100
                        addAmino2dict(Adomain,aminoName,score)

    # print per_orf_nrpspredictor2_scores_SVM
    per_orf_nrpspredictor2_stachsvm_aminos = {}
    per_orf_nrpspredictor2_stachsvm_score = {}
    per_orf_nrpspredictor2_stachsvm_rank = {}
    per_orf_nrpspredictor2_stachsvm_AAname = {}

    for orf in sorted(per_orf_nrpspredictor2_SVMs):
        addedaminos = []
        if orf in per_orf_nrpspredictor2_SVMs:
            if orf not in per_orf_nrpspredictor2_stachsvm_aminos:
                per_orf_nrpspredictor2_stachsvm_aminos[orf] = {}
            for i in range(len(per_orf_nrpspredictor2_SVMs[orf])):
                amino = per_orf_nrpspredictor2_SVMs[orf][i]
                aminoName = per_orf_nrpspredictor2_AAname_SVM[orf][i]

                if amino in addedaminos:
                    continue
                addedaminos.append(amino)
                # stachcode_score = per_orf_nrpspredictor2_scores[orf][0][i]
                if amino in per_orf_nrpspredictor2_SVMs[orf]:
                    j = per_orf_nrpspredictor2_SVMs[orf].index(amino)
                    svm_score = per_orf_nrpspredictor2_scores_SVM[orf][j]
                    new_score = svm_score
                    # new_score = (svm_score + stachcode_score)/2.0
                else:
                    new_score = (0 + stachcode_score)/2.0

                if new_score in per_orf_nrpspredictor2_stachsvm_aminos[orf]:
                    per_orf_nrpspredictor2_stachsvm_aminos[orf][new_score].append( (amino,aminoName) )
                else:
                    per_orf_nrpspredictor2_stachsvm_aminos[orf][new_score] = [(amino,aminoName)]

    per_ctg_nrpspred_byrank = {} # keys are contigs, the second key is the rank, third key is A-domain (orf), then it's a list of triples (amino,score,rank,aminoName)
    for orf in sorted(per_orf_nrpspredictor2_stachsvm_aminos):
        rank = -1 
        ctg = orf.split("_")[0]

        if ctg not in per_ctg_nrpspred_byrank:
            per_ctg_nrpspred_byrank[ctg] = {}
            for r in range(3):
                per_ctg_nrpspred_byrank[ctg][r] = {}

        if not per_orf_nrpspredictor2_stachsvm_aminos[orf]:
            continue
        maxscore = max(per_orf_nrpspredictor2_stachsvm_aminos[orf].keys())
        for score in sorted(set(per_orf_nrpspredictor2_stachsvm_aminos[orf]),reverse=True)[0:3]:
            # percent_max_score = score
            if score<50:
                continue
            percent_max_score = round(score/(maxscore*1.0),2)*100
            rank += 1
            if rank==3:
                break
            per_ctg_nrpspred_byrank[ctg][rank][orf] = []
            for (amino,aminoName) in per_orf_nrpspredictor2_stachsvm_aminos[orf][score]:
                per_ctg_nrpspred_byrank[ctg][rank][orf].append( (amino, percent_max_score, rank,aminoName) )

    final_nrps_info = dict()
    max_mass_threshold = 2000
    for possibleMass in [round(p/100.0,2) for p in range(30000, max_mass_threshold*100)]:
        final_nrps_info[round(possibleMass,2)] = {}

    
    for ctg in per_ctg_nrpspred_byrank:
        if ctg in per_ctgs_NRP_groups_domains:
            ctg_nrps_info = final_putative_NRP_gbk_code(
                per_ctg_nrpspred_byrank[ctg],
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
    def createPNPGraphLinesLinear(pnp,name):

        num = len(pnp)
        lines= []
        nameLine = name
        lines.append(nameLine) 
        newline = "number of components : "+ str(num)
        lines.append(newline)
        for i in range(num):
            aa = pnp[i]
            if i==0:
                aa= pnp[i]+1.007
            if i == num-1 :
                aa = pnp[i]+18.01528
            newline = str(i) + " CXHXNX " + str(aa)            
            lines.append(newline)
        newline = "number of bonds : " + str(num-1)
        lines.append(newline)
        for i in range(num-1):
            newline = str(i) + " -NC> " + str(i+1)
            lines.append(newline)
        lines.append(newline)
        return lines

    scores = set()
    
    peptidesUptoCorrectScore = []
    correctFound = False
    foundOne = False

    for nrpPair in candidateNRPSs:
        nrps = nrpPair[1]
        

        compoundName = nrpPair[0]

        graphFile = open(output+"_candidateNRPs.GRAPHS.txt","a")
        #'''

        linesToAdd = createPNPGraphLines(nrps,compoundName)
        foundOne = True
        for line in linesToAdd:
            graphFile.write(line+"\n")
        #'''linesToAdd = []
        origname = compoundName
        compoundName += "lin1"
        compoundName = origname + "_lin1"
        linesToAdd = createPNPGraphLinesLinear(nrps,compoundName)
        for line in linesToAdd:
            graphFile.write(line+"\n")
        compoundName = origname + "_revlin"
        revnrps = tuple(list(nrps))
        linesToAdd = createPNPGraphLinesLinear(revnrps,compoundName)
        for line in linesToAdd:
            graphFile.write(line+"\n")
    return foundOne
