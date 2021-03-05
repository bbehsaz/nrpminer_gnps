# from read_gbk_files import *
import os
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import re

class Reverse:
    """Iterator for looping over a sequence backwards."""
    def __init__(self, data):
        self.data = data
        self.index = len(data)

    def __iter__(self):
        return self

    def next(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.data[self.index]


def nrpspredictor_clusters():
    """ using the gbk predictios read in NRPSpredictor3 """
    small_clusters= set()
    small_clusters_dict = {}
    small_clusters.add(tuple(['gly','ala']))
    small_clusters.add(tuple(['val','leu','ile']))
    small_clusters.add(tuple(['ser']))
    small_clusters.add(tuple(['thr']))
    small_clusters.add(tuple(['phe','trp']))
    small_clusters.add(tuple(['tyr','bht']))
    small_clusters.add(tuple(['asp','asn']))
    small_clusters.add(tuple(['glu','gln']))
    small_clusters.add(tuple(['orn']))
    small_clusters.add(tuple(['arg']))
    small_clusters.add(tuple(['pro']))
    for smallCluster in small_clusters:
        for amino in smallCluster:
            small_clusters_dict[amino] = smallCluster            

    large_clusters= set()
    large_clusters.add(tuple(['gly','ala','val','leu','ile']))    
    large_clusters.add(tuple(['ser','thr'])) 
    large_clusters.add(tuple(['phe','trp','tyr']))    
    large_clusters.add(tuple(['asp','asn','glu','gln']))
    large_clusters.add(tuple(['orn','lys','arg']))
    large_clusters.add(tuple(['cys','no_call']))
    large_clusters.add(tuple(['Pro']))
    large_clusters_dict = {}
    for largeCluster in large_clusters:
        for amino in largeCluster:
            large_clusters_dict[amino] = largeCluster

    return small_clusters_dict, large_clusters_dict


    



def readSVMs(svmsFileAddress,topAA,minScore,cyclominerPath,antismash_res_dir):
    nrpspredictor_version = 3
    if "nrpspredictor3" in svmsFileAddress:
        nrpspredictor_version = 3
    if nrpspredictor_version == 2: 
        finalAntiSmashSVMPredicts,finalAntiSmashSVMPredictsScore,finalAntiSmashSVMPredictsRank,finalAntiSmashSVMPredictsAAname = nrpsredictor2_SVM(svmsFileAddress,topAA,minScore,cyclominerPath)
    else:
        finalAntiSmashSVMPredicts,finalAntiSmashSVMPredictsScore,finalAntiSmashSVMPredictsRank,finalAntiSmashSVMPredictsAAname = nrpsredictor2_SVM(svmsFileAddress,topAA,minScore,cyclominerPath)
    return finalAntiSmashSVMPredicts,finalAntiSmashSVMPredictsScore,finalAntiSmashSVMPredictsRank,finalAntiSmashSVMPredictsAAname



def nrpsredictor2_SVM(svmsFileAddress,topAA,minScore,cyclominerPath):
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashSVMPredicts ={}
    antiSmashSVMPredictsScore = {}
    antiSmashSVMPredictsRank = {}
    antiSmashSVMPredictsAAname = {}

    import numpy as np
    svmFile = open(svmsFileAddress,"r")
    line =svmFile.readline().strip()
    # new_Adomain_name = line.split()[0] #sequence name 
    while(True):
        line =svmFile.readline().strip()
        if not line: 
            break
        lineSplit = line.split()
        Adomain_name = lineSplit[0]
        predictions = [x.split(",") for x in line.split()[4:7]]

        
        if Adomain_name not in antiSmashSVMPredicts:
            antiSmashSVMPredicts[Adomain_name] = []
            antiSmashSVMPredictsScore[Adomain_name] = []
            antiSmashSVMPredictsRank[Adomain_name] = []
            antiSmashSVMPredictsAAname[Adomain_name] = []
        #score the SVM predcitions based on the cluster they appear in
        for aminoName in predictions[2]:
            if aminoName not in aminFullName2Mass:
                continue
            score = 100
            rank = 0
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        for aminoName in predictions[1]:
            if aminoName in predictions[2]:
                continue
            if aminoName not in aminFullName2Mass:
                continue
            score = 90
            rank = 1
            if len(antiSmashSVMPredicts)==0:
                score = 100
                rank = 0
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        for aminoName in predictions[0]:
            if aminoName in predictions[1]:
                continue
            if aminoName not in aminFullName2Mass:
                continue
            score = 80
            
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)

    finalAntiSmashSVMPredicts = {}
    finalAntiSmashSVMPredictsScore = {}
    finalAntiSmashSVMPredictsRank = {}
    finalAntiSmashSVMPredictsAAname = {}
    for contig in antiSmashSVMPredicts:
        if len(antiSmashSVMPredicts[contig]) <0:
            continue
        finalAntiSmashSVMPredicts[contig] = antiSmashSVMPredicts[contig]
        finalAntiSmashSVMPredictsScore[contig] = antiSmashSVMPredictsScore[contig]
        finalAntiSmashSVMPredictsRank[contig] = antiSmashSVMPredictsRank[contig]
        finalAntiSmashSVMPredictsAAname[contig] = antiSmashSVMPredictsAAname[contig]

    return finalAntiSmashSVMPredicts,finalAntiSmashSVMPredictsScore,finalAntiSmashSVMPredictsRank,finalAntiSmashSVMPredictsAAname

def readCodes(codesFileAddress,topAA,minScore,cyclominerPath): ##Read the NRPSpredictor2 code files. we do not filter at this point on any point. JUST READ the codes file.
    from database_functions import getBuildingBlocksForNRPSprediction
    import numpy as np
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashPredicts ={}
    antiSmashPredictsScore = {}
    antiSmashPredictsRank = {}
    antiSmashPredictsAAname = {}
    with open(codesFileAddress,"r") as gbkFile:
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
            antiSmashPredictsAAname[contig] = {}
            while new_contig ==contig:
                antiSmashPredicts[contig][pos] = []
                antiSmashPredictsScore[contig][pos] = []
                antiSmashPredictsRank[contig][pos] = []
                antiSmashPredictsAAname[contig][pos] = []
                old_score =110
                rank=-1

                for i in range(30):
                    # if len(set(antiSmashPredictsScore[contig][pos]))>topAA:
                    #     break
                    aminoScore = aminoAcids[i]
                    #here is where the amino acids are read
                    aminoName = aminoScore.split("(")[0]

                    score = float(aminoScore.split("(")[1].split(")")[0])
                    if score<old_score:
                            # if len(antiSmashPredicts[contig][pos])>=topAA:
                            #     break
                            rank+=1
                            if len(set(antiSmashPredictsScore[contig][pos]))>=topAA:
                                break
                            if rank==3:
                                break
                    if i==0:
                        adomainmaxscore = score
                        #read the amino acids predicted in the code file
                    if aminoName in aminFullName2Mass:
                        if score<minScore:
                            continue
                        antiSmashPredictsAAname[contig][pos].append(aminoName)
                        # print "hello?"
                        # print aminoName
                        antiSmashPredicts[contig][pos].append(aminFullName2Mass[aminoName])
                                                
                        newscore = score

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
                new_contig =  lineSplit[0]
                predictions = line.split()[2]
                aminoAcids = predictions.split(";")

    finalAntiSmashPredicts = {}
    finalAntiSmashPredictsScore = {}
    finalAntiSmashPredictsRank = {}
    finalAntiSmashPredictsAAname = {}
    for contig in antiSmashPredicts:
        if len(antiSmashPredicts[contig]) <0:
            continue
        finalAntiSmashPredicts[contig] = antiSmashPredicts[contig]
        finalAntiSmashPredictsScore[contig] = antiSmashPredictsScore[contig]
        finalAntiSmashPredictsRank[contig] = antiSmashPredictsRank[contig]
        finalAntiSmashPredictsAAname[contig] = antiSmashPredictsAAname[contig]

    
    return finalAntiSmashPredicts,finalAntiSmashPredictsScore,finalAntiSmashPredictsRank, finalAntiSmashPredictsAAname

#read the NRPS monomers
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



def nrpsredictor3_SVM(svmsFileAddress,topAA,minScore,cyclominerPath):
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashSVMPredicts ={}
    antiSmashSVMPredictsScore = {}
    antiSmashSVMPredictsRank = {}
    antiSmashSVMPredictsAAname = {}

    import numpy as np

    svmFile = open(svmsFileAddress,"r")
    line =svmFile.readline().strip()
    # new_Adomain_name = line.split()[0] #sequence name 
    while(True):
        line =svmFile.readline().strip()
        if not line: 
            break
        if "single" in line:
            continue
        lineSplit = line.split()
        Adomain_name = "_".join(lineSplit[0].split("_")[0:3])

        predictions = [x.split(",") for x in line.split()[4:7]]  

        antiSmashSVMPredicts[Adomain_name] = []
        antiSmashSVMPredictsScore[Adomain_name] = []
        antiSmashSVMPredictsRank[Adomain_name] = []
        antiSmashSVMPredictsAAname[Adomain_name] = []

        # print predictions
        #score the SVM predcitions based on the cluster they appear in
        if len(predictions) >2:
            for aminoName in predictions[2]:
                if aminoName in antiSmashSVMPredictsAAname[Adomain_name]:
                    continue
                if aminoName not in aminFullName2Mass:
                    continue
                score = 100
                rank = 0
                antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
                antiSmashSVMPredictsScore[Adomain_name].append(score)
                antiSmashSVMPredictsRank[Adomain_name].append(rank)
                antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        if len(predictions)>1:
            for aminoName in predictions[1]:
                if aminoName in antiSmashSVMPredictsAAname[Adomain_name]:
                    continue
                # if len(predictions)>2:
                #     if aminoName in predictions[2]:
                #         continue
                if aminoName not in aminFullName2Mass:
                    continue
                score = 90
                rank = 1
                if len(antiSmashSVMPredicts)==0:
                    score = 100
                    rank = 0
                antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
                antiSmashSVMPredictsScore[Adomain_name].append(score)
                rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
                antiSmashSVMPredictsRank[Adomain_name].append(rank)
                antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        if len(predictions)>0:
            for aminoName in predictions[0]:
                if aminoName in antiSmashSVMPredictsAAname[Adomain_name]:
                    continue
                # if len(predictions)>1:
                #     if aminoName in predictions[1]:
                #         continue
                if aminoName not in aminFullName2Mass:
                    continue
                score = 80
                
                antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
                antiSmashSVMPredictsScore[Adomain_name].append(score)
                rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
                antiSmashSVMPredictsRank[Adomain_name].append(rank)
                antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        ''' uncomment for the antismash nrpspredictor2
        predictions = [x.split(",") for x in line.split()[8:11]]        
        # antiSmashSVMPredicts, antiSmashSVMPredictsScore, antiSmashSVMPredictsRank, antiSmashSVMPredictsAAname = read_SVM_predictions_line_withScores(antiSmashSVMPredicts, predictions, aminFullName2Mass)
        #     if Adomain_name not in antiSmashSVMPredicts:

        print predictions
        #score the SVM predcitions based on the cluster they appear in
        if len(predictions) == 0 :
            continue
        for aminoName in predictions[2]:
            if aminoName not in aminFullName2Mass:
                continue
            score = 100
            rank = 0
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        for aminoName in predictions[1]:
            if aminoName in predictions[2]:
                continue
            if aminoName not in aminFullName2Mass:
                continue
            score = 90
            rank = 1
            if len(antiSmashSVMPredicts)==0:
                score = 100
                rank = 0
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        for aminoName in predictions[0]:
            if aminoName in predictions[1]:
                continue
            if aminoName not in aminFullName2Mass:
                continue
            score = 80
            
            antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashSVMPredictsScore[Adomain_name].append(score)
            rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
            antiSmashSVMPredictsRank[Adomain_name].append(rank)
            antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
        ''' #uncomment for the antismash nrpspredictor2
    # print antiSmashSVMPredicts
    # print antiSmashSVMPredictsScore
    finalAntiSmashSVMPredicts = {}
    finalAntiSmashSVMPredictsScore = {}
    finalAntiSmashSVMPredictsRank = {}
    finalAntiSmashSVMPredictsAAname = {}
    for contig in antiSmashSVMPredicts:
        if len(antiSmashSVMPredicts[contig]) <0:
            continue
        finalAntiSmashSVMPredicts[contig] = antiSmashSVMPredicts[contig]
        finalAntiSmashSVMPredictsScore[contig] = antiSmashSVMPredictsScore[contig]
        finalAntiSmashSVMPredictsRank[contig] = antiSmashSVMPredictsRank[contig]
        finalAntiSmashSVMPredictsAAname[contig] = antiSmashSVMPredictsAAname[contig]
    return finalAntiSmashSVMPredicts,finalAntiSmashSVMPredictsScore,finalAntiSmashSVMPredictsRank,finalAntiSmashSVMPredictsAAname


def read_SVM_predictions_line_withScores(antiSmashSVMPredicts, predictions, aminFullName2Mass):
    if Adomain_name not in antiSmashSVMPredicts:
        antiSmashSVMPredicts[Adomain_name] = []
        antiSmashSVMPredictsScore[Adomain_name] = []
        antiSmashSVMPredictsRank[Adomain_name] = []
        antiSmashSVMPredictsAAname[Adomain_name] = []
    #score the SVM predcitions based on the cluster they appear in
    for aminoName in predictions[2]:
        if aminoName not in aminFullName2Mass:
            continue
        score = 100
        rank = 0
        antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
        antiSmashSVMPredictsScore[Adomain_name].append(score)
        antiSmashSVMPredictsRank[Adomain_name].append(rank)
        antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
    for aminoName in predictions[1]:
        if aminoName in predictions[2]:
            continue
        if aminoName not in aminFullName2Mass:
            continue
        score = 90
        rank = 1
        if len(antiSmashSVMPredicts)==0:
            score = 100
            rank = 0
        antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
        antiSmashSVMPredictsScore[Adomain_name].append(score)
        rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
        antiSmashSVMPredictsRank[Adomain_name].append(rank)
        antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
    for aminoName in predictions[0]:
        if aminoName in predictions[1]:
            continue
        if aminoName not in aminFullName2Mass:
            continue
        score = 80
        
        antiSmashSVMPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
        antiSmashSVMPredictsScore[Adomain_name].append(score)
        rank = len(set(antiSmashSVMPredictsScore[Adomain_name]))
        antiSmashSVMPredictsRank[Adomain_name].append(rank)
        antiSmashSVMPredictsAAname[Adomain_name].append(aminoName)
    return antiSmashSVMPredicts, antiSmashSVMPredictsScore, antiSmashSVMPredictsRank, antiSmashSVMPredictsAAname


def readINDs(indFileAddress,topAA,minScore,cyclominerPath):
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashINDPredicts = {}
    antiSmashINDPredictsScore = {}
    antiSmashINDPredictsRank = {}
    antiSmashINDPredictsAAname = {}

    import numpy as np
    indFile = open(indFileAddress,"r")
    line =indFile.readline().strip()
    # new_Adomain_name = line.split()[0] #sequence name 
    while(True):
        line =indFile.readline().strip()        
        if not line: 
            break
        lineSplit = line.split()
        Adomain_name = "_".join(lineSplit[0].split("_")[0:3])
        if Adomain_name not in antiSmashINDPredicts:
            antiSmashINDPredicts[Adomain_name] = []
            antiSmashINDPredictsScore[Adomain_name] = []
            antiSmashINDPredictsRank[Adomain_name] = []
            antiSmashINDPredictsAAname[Adomain_name] = []
        aminoName = lineSplit[2]
        if "ASM" in line or "SVM" in line or "pHMM" in line:
            if aminoName not in aminFullName2Mass:
                continue
            if aminoName in antiSmashINDPredictsAAname[Adomain_name]:
                continue
            score = 100
            rank = 0
            antiSmashINDPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashINDPredictsScore[Adomain_name].append(score)
            antiSmashINDPredictsRank[Adomain_name].append(rank)
            antiSmashINDPredictsAAname[Adomain_name].append(aminoName)

    # print "ola?"
    # print antiSmashINDPredicts
    # print antiSmashINDPredictsScore
    finalantiSmashINDPredicts = {}
    finalantiSmashINDPredictsScore = {}
    finalantiSmashINDPredictsRank = {}
    finalantiSmashINDPredictsAAname = {}
    for contig in antiSmashINDPredicts:
        if len(antiSmashINDPredicts[contig]) <0:
            continue
        finalantiSmashINDPredicts[contig] = antiSmashINDPredicts[contig]
        finalantiSmashINDPredictsScore[contig] = antiSmashINDPredictsScore[contig]
        finalantiSmashINDPredictsRank[contig] = antiSmashINDPredictsRank[contig]
        finalantiSmashINDPredictsAAname[contig] = antiSmashINDPredictsAAname[contig]

    return finalantiSmashINDPredicts,finalantiSmashINDPredictsScore,finalantiSmashINDPredictsRank,finalantiSmashINDPredictsAAname




def get_specificty_predictions_genbank_features(
    gb_record,
    feature_type,
    qualifier,
    predictor) :
    # fiven a gbk record finds the specificity prediction based on predictor 
    # predictor between ['Stachelhaus','NRPSpredictor3','pHMM','SANDPUMA'])
    import os
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    import re
    answer = 'no_call'
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :

                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for prediction in feature.qualifiers[qualifier]:
                    if predictor in prediction:
                        # print prediction
                        # print predictor
                        # print feature_type
                        # print qualifier

                        answer = prediction.split(":")[1][1:]
                # for value in feature.qualifiers[qualifier] :
                #     if value in answer :
                #         print "WARNING - Duplicate key %s for %s features %i and %i" \
                #            % (value, feature_type, answer[value], index)
                #     else :
                #         answer[value] = index
    return answer.split("|")



def readGBK_predictions(gb_record,topAA,minScore,cyclominerPath):
    from database_functions import get_specificty_predictions_genbank_features
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashGBKPredicts = {}
    antiSmashGBKPredictsScore = {}
    antiSmashGBKPredictsRank = {}
    antiSmashGBKPredictsAAname = {}
    
    import numpy as np
    indFile = open(indFileAddress,"r")
    line =indFile.readline().strip()
    # new_Adomain_name = line.split()[0] #sequence name 

    while(True):
        line =indFile.readline().strip()        
        if not line: 
            break
        lineSplit = line.split()
        Adomain_name = "_".join(lineSplit[0].split("_")[0:3])
        if Adomain_name not in antiSmashINDPredicts:
            antiSmashINDPredicts[Adomain_name] = []
            antiSmashINDPredictsScore[Adomain_name] = []
            antiSmashINDPredictsRank[Adomain_name] = []
            antiSmashINDPredictsAAname[Adomain_name] = []
        aminoName = lineSplit[2]
        if "ASM" in line or "SVM" in line or "pHMM" in line:
            if aminoName not in aminFullName2Mass:
                continue
            if aminoName in antiSmashINDPredictsAAname[Adomain_name]:
                continue
            score = 100
            rank = 0
            antiSmashINDPredicts[Adomain_name].append(aminFullName2Mass[aminoName])
            antiSmashINDPredictsScore[Adomain_name].append(score)
            antiSmashINDPredictsRank[Adomain_name].append(rank)
            antiSmashINDPredictsAAname[Adomain_name].append(aminoName)

    finalantiSmashINDPredicts = {}
    finalantiSmashINDPredictsScore = {}
    finalantiSmashINDPredictsRank = {}
    finalantiSmashINDPredictsAAname = {}
    for contig in antiSmashINDPredicts:
        if len(antiSmashINDPredicts[contig]) <0:
            continue
        finalantiSmashINDPredicts[contig] = antiSmashINDPredicts[contig]
        finalantiSmashINDPredictsScore[contig] = antiSmashINDPredictsScore[contig]
        finalantiSmashINDPredictsRank[contig] = antiSmashINDPredictsRank[contig]
        finalantiSmashINDPredictsAAname[contig] = antiSmashINDPredictsAAname[contig]

    return finalantiSmashINDPredicts,finalantiSmashINDPredictsScore,finalantiSmashINDPredictsRank,finalantiSmashINDPredictsAAname








