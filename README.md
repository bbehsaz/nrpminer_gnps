# NRPminer 1.0: Scalable Peptidogenomics Algorithm for Non-Ribosomal Peptide Discovery 

NRPminer is a modification-tolerant scalable metabologenomics method for discovering NRPs by integrating (meta)genomics and tandem mass spectrometry (metabolomics) data

## Dependencies (add binaries to PATH):

	- NPDtools (https://github.com/ablab/npdtools)

## System requirements:

	- Linux or macOS
	- Python 2.7
	- GNU sed 

## To configure/test run:

     make Makefile

## Usage examples: 


	nrpminer.py -s [spectrum_file]  --antismash_resgbk [path_to_antiSMASH_generated_gbk] --orfDel [num_orfs_deleted] --derepdir [path_to_NPDtools_bin] --maxmod [maximum_modification_mass] -o [output_folder] --pvalue [pvalue_threshold] --threads [number of threads] --nrpspks_predictions_txt [path to SVM and CODE predictions folder]       



For the full list of available options please run

     nrpminer.py -h


Output:

** all_significant_psms.NRPminer.results.tsv ::            list of all PSMs with significant PSM P-value as well as the sequence of the identified NRP

