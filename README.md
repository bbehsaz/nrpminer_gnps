# NRPminer 0.2: Scalable Peptidogenomics Algorithm for Non-Ribosomal Peptide Discovery 

CycloNovo is a new algorithm that identifies spectra generated from cyclopeptides in large mass spectrometry datasets. CycloNovo can also de novo sequence the cyclopeptides represented by identified cyclospectra.


Dependencies (add binaries to PATH):

	- NPDtools (https://github.com/ablab/npdtools)

System requirements:

	- Linux or macOS
	- Python 2.7
	- GNU sed 

To configure/test run:

     make Makefile

Usage examples: 


	nrpminer.py -s [spectrum_file]  --antismash_resdir [path_to_antiSMASH_results_folder] --orfDel [num_orfs_deleted] --derepdir [path_to_NPDtools_bin] --maxmod [maximum_modification_mass] -o [output_folder] --pval_threshold [pvalue_threshold]      


In new version, NRPminer only requires gbk files generated by antiSMASH (in 'path_to_antiSMASH_results_folder').

For the full list of available options please run

     nrpminer.py -h


Output:

* nrpminer_significant_psms.txt:                    list of all PSMs with significant PSM P-value as well as the sequence of the identified NRP
* nrpminer_all_psms.txt:                    list of all PSMs as well as the sequence of the identified NRP

