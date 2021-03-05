

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
def readSVMs(svmsFileAddress,topAA,minScore,cyclominerPath):
    standardMasses,aminFullName2Mass = getBuildingBlocksForNRPSprediction(cyclominerPath+"/configs/nrp_hosein_mass.txt")
    antiSmashSVMPredicts ={}
    antiSmashSVMPredictsScore = {}
    antiSmashSVMPredictsRank = {}
    import numpy as np
    svmFile = open(svmsFileAddress,"r")
    line =svmFile.readline().strip()
#     new_contig = line.split()[0].split("_A")[0] #these are orfs not contigs
    new_contig = line.split()[0] #sequence name 
    svmpreds = [line.split("\t") for line in svmFile.read().split("\n") if "\t" in line]
    for prediction in svmpreds:
        NRPSname = prediction[0].rpartition("_")[0]
        if not nrps_info[idx]['svmpredictions'].has_key(NRPSname):
            nrps_info[idx]['svmpredictions'][NRPSname] = {}
        domainname = prediction[0].rpartition("_")[2]
        nrps_info[idx]['svmpredictions'][NRPSname][domainname] = prediction[3:7]

