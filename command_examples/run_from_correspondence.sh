#######################
# Script that automatically runs the entire datasets using their correspondence files. 
# Inputs:
#(1) correspondence_file.csv 
#(2) folder containing all the MS/MS files 
#(3) folder containing the antismash results 
#(4) projectname (will be used as name of outputfolder to store the entire run) 
#(5) a path to the folder contaning NPDtools/VarQuest executables
#######################

xaa=$1 
ms_folder=$2
antismash_folder=$3
outputfolder=$4
dereppath=$5
nrpminerpath=$6

mkdir -p $outputfolder

python "$nrpminerpath"/scripts/process_correspondence.py "$xaa" $outputfolder

if [ -f "$outputfolder"/NRPminer_formatted_correspondence.txt ]; then
	while IFS=$'\t' read -r -a myArray
		do
			genome=$(echo "${myArray}")
			genomefile=$(find "$antismash_folder" -name "*$genome*" -type d -maxdepth 1)
			echo "==============================="
			echo "Processing ... $genomefile"
			python "$nrpminerpath"/nrpminer.py --delete_orfs 2 --derepdir "$dereppath" --maxmod 150 --topAA 2 -o "$outputfolder"/NRPDBs/"$genome" --delete_aas 0 --makeNRPDBonly --antismash_resdir "$genomefile" > "$outputfolder"/"$genome"_orfdel0.log.txt
			wait 
			for specname in "${myArray[@]:0}"
				do
					find "$ms_folder" -name "*$specname*" -type f -maxdepth 2 | while read specfile;
					do
						find "$outputfolder"/NRPDBs/"$genome" -type d -name "NRPminer_graphs" | while read NRP_DB;
						do 
							#python "$nrpminerpath"/nrpminer.py -s "$specfile" --delete_orfs 0 --derepdir "$dereppath" --maxmod 150 --topAA 2 -o "$outputfolder"/spectrum_nrp_matches/"$specname"_vs_"$genome" --coreNRPDB "$NRP_DB" --coreNRPDBlib "$NRP_DB"/library.info.graphs 
							wait 
						done
						wait
					done
					wait
				done
				wait
		done < "$outputfolder"/NRPminer_formatted_correspondence.txt
	wait
fi


#for i in {0..2}
#do
# # your-unix-command-here
# echo $i
# python nrpquest_nrpmasses_allderep.py -s $mgf --antismash_resdir $antismashres --delete_orfs $i --derepdir ./../../dereplicator/ --maxmod 150 --topAA 2 -o test/"$outname"/orfs_deleted_"$i" --delete_aas 0
#done 

