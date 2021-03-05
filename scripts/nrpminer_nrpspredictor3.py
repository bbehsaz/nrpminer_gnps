#!/usr/bin/env python
#################################################
# NRPminer: algorithm for high-throughput integration of metabolomics and (meta)genomics for non-ribosomal peptide discovery
# Developed by Bahar Behsaz (June 2018)
# Pevzner lab 
# University of California at San Diego, La Jolla, CA, USA
#################################################
import sys
import os
import signal

import traceback
import errno
import logging
from os.path import join, abspath, isdir
from datetime import datetime
from site import addsitedir
from optparse import OptionParser
import argparse
import os


from scripts.readMGF import readMGF
from scripts.readStream import readStream
from scripts.read_string_mgf_vectored import readStringMGF

from scripts.create_putative_nrp_database_nrpspredictor3 import *

from scripts.read_nrpspredictor2_codes import readCodes
from scripts.read_nrpspredictor2_codes import getBuildingBlocksForNRPSprediction

from scripts.spectral_convolutions_analysis import generate_aa_convolutions_vectorized
from scripts.spectral_convolutions_analysis import convolution
from scripts.createNRPdatabase_contig_score_based import generateNRPS3mers
from scripts.createNRPdatabase_contig_score_based import generateNRPSpredict
from scripts.createNRPdatabase_contig_score_based import output_nrpsgraph_bcyclic
from scripts.create_putative_nrp_database import output_nrpsgraph
from scripts.createNRPdatabase_contig_score_based import findNRPPredicts

from utils.sharedFunctions import getBuildingBlocks
from utils.sharedFunctions import writeOriginalSpectra
from utils.sharedFunctions import findRealPepMass



from utils.spectra_utils import get_spectra_fpaths
from utils.file_utils import verify_file
from utils.file_utils import is_valid_file
from utils.file_utils import CleanChildProcesses
from utils.file_utils import readable_dir

from utils.common import parse_params_xml
from utils.file_utils import removeDir
from utils.file_utils import mkdir_p


import math
from operator import itemgetter
import itertools
from subprocess import Popen,PIPE
 
import Bio
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import shutil

def get_parser():
    """Parse arguments and check if the spectrum file exists. """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(prog="NRPminer.py",description= "Algorith for analyzing cyclopeptides using Tandem Mass Spectra ...", 
        formatter_class=ArgumentDefaultsHelpFormatter)

    ### metabolomics dataset to search against the database of putative NRPs 
    parser.add_argument("-s", "--spectrum", dest="mgfFile", type=lambda x: is_valid_file(parser, x), help="Input spectra (mgf/mzXML)", metavar="FILE",default="")

    # parser.add_argument("--nrpspredictor_code", dest="nrpspredictor_code", type=lambda x: is_valid_file(parser, x), help="Input nrpspredictor result (*codes.txt) ", metavar="FILE", default='null')
    parser.add_argument('--antismash_resdir', help="antiSMASH result diroctory ", default="")
    parser.add_argument('--nrpspks_predictions_txt', help="Directory including the Stachelhaus and SVM prediction TXT files. Previously nrpspks_predictions_txt folder. ", default="")
    parser.add_argument("--threads", dest='numthreads', help='Number of threads', default=int(1),type=int)

    ### make a putative database of putative NRPs 
    parser.add_argument("-o", "--out_dir", dest="outdir", help="Output dir",required=True)
    parser.add_argument("--pname", dest="pname", help="Project name used as prefix for generated files", default="NRPminer")
    parser.add_argument("--derepdir", dest="derepdir", help="Directory containing Dereplicator master ", required=True)
    parser.add_argument("--maxmod", dest="maxmod", help='Maximum mass modification to consider', default=150,type=float)
    parser.add_argument("--minAdomainscore", dest="minAdomainscore", help='NRPs with minimum A-domain score', default=50,type=int)
    parser.add_argument("--topAA", dest="topAA", help='Maximum rank of amino acids considerd for each A-domain', default=2,type=int)
    parser.add_argument("--orfDel", dest="delete_orfs", help='Number of orfs deleted', default=0,type=int)
    parser.add_argument("--orfDup", dest="duplicate_orfs", help='Number of times an ORF is duplicated', default=0,type=int)
    parser.add_argument("--N", dest="thresholdN", help='Number of orfs deleted', default=1000,type=int)
    parser.add_argument("--delete_aas", dest="delete_aminos", help='Number of amino acids deleted', default=0,type=int)
    parser.add_argument("--aa_rankThreshold", dest="aa_rankThreshold", help='Maximum sum of ranks in the putative NRPs', default=2,type=int)
    parser.add_argument('--peptideonly', dest='peponly', help="Peptide spectra only", action='store_true')
    parser.add_argument('--makeNRPDBonly', dest='makeNRPDBonly', help="Peptide spectra only", action='store_true',default=False)

    ### use pre made putative NRP database
    parser.add_argument('--coreNRPDB', dest="coreNRPDB", help="Absolute path to preprocessed putative NRP database (one before graph_dir)", default="")
    parser.add_argument('--coreNRPDBlib', dest="coreNRPDBlib", help="Absolute path to the library of preprocessed putative NRPs in the input database (library.info.graphs)", default="")
    
    parser.add_argument("-d", "--precursor_ion_thresh", dest='delta', help='Precursor Ion Thresh for custom Running mode', default='0.02')
    parser.add_argument("-e","--fragment_ion_thresh", dest='e', help='Fragment Ion Thresh for custom Running mode', default='0.02')
    parser.add_argument("-v","--verbosity", dest="verbose", action="store_true")

    #p-value threshold 
    parser.add_argument("--pvalue", dest="pvalue", help='Maximum mass modification to consider', default=1.0e-15,type=float)
    return parser


def get_params_and_mapping(opts):
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  # TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping

if __name__ == "__main__":
    args = get_parser().parse_args()
    NRPminerPath = os.path.dirname(os.path.realpath(__file__))
    if args.verbose:
        def verboseprint(stringtoprint):
            print(stringtoprint)
    else:
        verboseprint = lambda *a, **k: None # do-nothing function
    verboseprint("####################################################################")
    verboseprint("++ Running NRPminer from directory: " + NRPminerPath)
    if len(args.antismash_resdir) ==0 and len(args.coreNRPDB)==0:
        print "Need either an antiSMASH result directory or a core NRP database generated by NRPminer previously"
        exit()
    elif len(args.antismash_resdir) > 0 and len(args.coreNRPDB) > 0:
        print "Need either an antiSMASH result directory OR a core NRP database generated by NRPminer previously; but not both."
        exit()

    if (not args.makeNRPDBonly) and len(args.mgfFile)==0:
        print "No input specturm (metabolites file). "
        exit()
    if args.makeNRPDBonly and len(args.coreNRPDB)>0:
        print "Already input NRP DB, not making NRP DB"
        exit()

    verboseprint("++ Reading input spectrum: " + args.mgfFile)
    verboseprint("++ Reading antismash predictions from dir: " + args.antismash_resdir)
    protonMass = 1.00728 #This is mass of hydrogen --- do not remove!
    digits = 3
    print args.pvalue
    exit()
    antiSMASH_directory_name =os.path.dirname(args.antismash_resdir)
    building_blocks_main = getBuildingBlocks(NRPminerPath+"/configs/nrp_hosein_mass.txt") #Dic all possible NRP amino acids
    building_blocks_stnd = getBuildingBlocks(NRPminerPath+"/configs/aminoacidMasses.txt") #Dic all possible stnd aa's amino acids
    polymer_repeat_units = getBuildingBlocks(NRPminerPath+"/configs/polymer_repeat_masses.txt") # Dic all polymer repeat units    
    shutil.rmtree(args.outdir, ignore_errors=True) #Delete any directory locatd at output directory
    
    #creating the output directory
    mkdir_p(args.outdir)
    mkdir_p(args.outdir+"/genomic_work")
    output = args.outdir + "/" + args.pname
    outtemp = args.outdir +"/genomic_work"+ "/" + "NRPminer"

    #check the input based on the chosen options:
    if args.peponly:
        peaksnIntensity,pepMasses,charges,retentions,fileLines, allSpectraVector_dic= readStringMGF(peaksfile)
    else:
        peaksnIntensity = {}
    if len(args.coreNRPDB)>0:
        if len(args.coreNRPDBlib)==0:
            print "library info for preprocessed core NRP database is required"
            exit()

    e = round(float(args.e),2)

    allFilteredSpecMasses = {}
    max_mass_threshold= 2000 #this is the maximum peptide mass accepted by NRPminer

    for possibleMass in [round(p/100.0,2) for p in range(30000, max_mass_threshold*100+1)]:
        if args.peponly:
            allFilteredSpecMasses[round(possibleMass,2)] = 0
        else:
            allFilteredSpecMasses[round(possibleMass,2)] = 1
    maxMod = args.maxmod

    
    rankThreshold = int(args.aa_rankThreshold)
    
    nameofcyclospecfile=str(args.mgfFile)
    cyclospectraAbsPath= os.path.abspath(nameofcyclospecfile)    
    dbpremade =False
    #### Use the already generated Core NRP database generated by NRPminer 
    if len(args.coreNRPDB)>0:   
        dbpremade = True
        nameofcoreNRPDB=str(args.coreNRPDB)
        nameofcoreNRPDBAbsPath= os.path.abspath(nameofcoreNRPDB)
        nameofcoreNRPDBLIB=str(args.coreNRPDBlib)
        nameofcoreNRPDBLIBAbsPath= os.path.abspath(nameofcoreNRPDBLIB)
        groupname = os.path.dirname(args.coreNRPDB+"/").split("/")[-3]
        mgffile  = os.path.dirname(args.mgfFile+"/").split("/")[-1]
        verboseprint("Using the following genome " + str(groupname))

    #### Use antiSMASH generated results to generated core NRPs ... 
    if not dbpremade:

        totalnum = 0
        a = 0
        groupname  = os.path.dirname(args.antismash_resdir+"/").split("/")[-1]
        if len(groupname)==0:
            print "The genome file has a problem"
            exit()

        outtemp_group = args.outdir +"/" + "genomic_work"+ "/" + groupname + "/" + args.pname         
        allnrpMasses = []
        numSpecMassMatched = 0
        linkedMasses = []
        thisbuildingblock = building_blocks_main.copy()
        final_numCandidatesPerSpectrum = {}
        verboseprint("Generating putative NRP database ...")
        final_nrps_info= create_putative_nrp_dataset(args.antismash_resdir, args.nrpspks_predictions_txt, allFilteredSpecMasses,NRPminerPath,rankThreshold,args.delete_orfs,args.delete_aminos,args.topAA,args.minAdomainscore,maxMod)
        nrps2check = []
        for mass in final_nrps_info:
                nrps2check += [(nrps[0],nrps[1]) for nrps in set(final_nrps_info[mass].items())]
        print "# of NRPminer generated putative NRPs: {0}".format(len(nrps2check))
        if len(nrps2check) < 1:
            print "No proper putative NRP found ... "
            exit(2)
        mkdir_p(args.outdir +"/" + "genomic_work"+ "/" + groupname)
        verboseprint("Writing the database ...")
        addMass= output_nrpsgraph(nrps2check,outtemp_group)
        totalnum = 0
        a = 0
        
    ##### sending spectra and the putative NRPs database to the modification tolerant search 
    nameofrecontfile = outtemp+"_spec_nrps_link.txt"
    massFileAbsPath = os.path.abspath(outtemp+"_masses.txt")    
    reconstFileAbsPath= os.path.abspath(nameofrecontfile)
    linkedMasses = [1]
    outtemp_group = args.outdir +"/" + "genomic_work"+ "/" + groupname + "/" + args.pname 
    nameofcyclospecfile = outtemp_group+"_cycloquest_nrps_spectra_linked.mgf"
    oldcyclospectraAbsPath= os.path.abspath(nameofcyclospecfile)
    outputDirAbsPath = os.path.dirname(oldcyclospectraAbsPath)
    outputWithAbsPath = outtemp_group
    derepDirAbsPath = os.path.dirname(os.path.abspath(args.derepdir+"/dereplicator.py"))
    finished = False
    permoutdir = os.path.dirname(output+"_NRPminer_summary.txt")

    with CleanChildProcesses():
        print args.makeNRPDBonly
        if args.makeNRPDBonly:
            proc = Popen(["sh " +NRPminerPath+'/scripts/run_format_DB.sh '+outputWithAbsPath+ " " + cyclospectraAbsPath + " " + str(NRPminerPath)  + " " + str(outputDirAbsPath) + " " + str(derepDirAbsPath) + " " + str(args.e)+ " " + str(int(args.maxmod)) + " " + str(permoutdir) + " "+ groupname + " " + str(args.numthreads)
            ] ,shell=True)
            try:
                proc.wait()
            except KeyboardInterrupt:
                try:
                    proc.terminate()
                except OSError:
                    pass
                    proc.wait

        elif len(args.coreNRPDB)>0:
            proc = Popen(["sh " +NRPminerPath+'/scripts/run_varquest_mass_nrp_onpremadeDB.sh '+outputWithAbsPath+ " " + cyclospectraAbsPath + " " + str(NRPminerPath)  + " " + str(outputDirAbsPath) + " " + str(derepDirAbsPath) + " " + str(args.e)+ " " + str(int(args.maxmod)) + " " + str(permoutdir) + " "+ groupname + " " + str(args.numthreads) + " " +str(nameofcoreNRPDBAbsPath) + " " +str(nameofcoreNRPDBLIBAbsPath) 
            ] ,shell=True)
            try:
                proc.wait()
            except KeyboardInterrupt:
                try:
                    proc.terminate()
                except OSError:
                    pass
                    proc.wait            
        else:

            proc = Popen(["sh " +NRPminerPath+'/scripts/run_varquest_mass_nrp.sh '+outputWithAbsPath+ " " + cyclospectraAbsPath + " " + str(NRPminerPath)  + " " + str(outputDirAbsPath) + " " + str(derepDirAbsPath) + " " + str(args.e)+ " " + str(int(args.maxmod)) + " " + str(permoutdir) + " "+ groupname + " " + str(args.numthreads)
            ] ,shell=True)
            try:
                proc.wait()
            except KeyboardInterrupt:
                try:
                    proc.terminate()
                except OSError:
                    pass
                    proc.wait





print "FINISHED"
exit()












