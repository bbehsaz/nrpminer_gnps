
import sys
import string 

def fix_psms(filepath):
	allnewlines = []
	with open(filepath.strip()) as psms_file:
		headerline = ""
		for line in psms_file:
			if "SpecFile" in line: 
				headerline = line.strip().split() 
				prefixhead = "\t".join(headerline[:5])
				newcolhead = "\t" + "ctgID+pos" + "\t"  + "nrpSeq" + "\t" + "nrpScore" + "\t"
				suffixhead = "\t".join(headerline[6:])
				newheaderline = prefixhead + newcolhead + suffixhead
				continue

			lineSplit = line.strip().split() 
			coltofix = lineSplit[5].split("_")
			nrpseq = coltofix[-1]
			nrplen = len(nrpseq.split("-"))
			# nrpNum = coltofix[-2]
			nrpscore = float(coltofix[-2])/(nrplen*1.0)
			# ctgseq = "_".join(coltofix[0:-3])
			prefixnewline = "\t".join(lineSplit[:5])

			ctgID = coltofix[0]
			minpos = coltofix[1]
			maxpos = coltofix[2]
			if coltofix[3] == "pos":				
				strand = "+"
			else:
				strand = "-"

			ctgseq = ctgID + "" + "("+minpos + ","+maxpos + "," + strand + ")"
			nrpseq_report = "cyc(" + nrpseq + ")"
			nrpseq_report = nrpseq
			newcols = "\t" + ctgseq + "\t" + nrpseq_report + "\t" + str(round(nrpscore,2)) + "\t"
			suffixnewline = "\t".join(lineSplit[6:])
			newline =prefixnewline + newcols + suffixnewline
			allnewlines.append(newline)
	newpsmsfilename = filepath.replace('.txt','') + ".NRPminer.results.tsv"
	with open(newpsmsfilename, "w") as fixedpsmsfile:
		fixedpsmsfile.write(newheaderline+"\n")
		for newline in allnewlines:
			fixedpsmsfile.write(newline+"\n")




fix_psms(sys.argv[1])
# fix_psms("/Users/bahar/workspace/npd_tools/nrpminer/nrpminer_antismash5/test_data/test_output_antismash5_2/psms_GCA_000719925.1_ASM71992v1_genomic_mod_all_significant_psms.txt")



