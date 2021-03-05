xaa=$1 
ms_folder=$2
antismash_folder=$3
outputfolder=$4
dereppath=$5
nrpminerpath=$6
orfdel=$7
#dereppath=/Users/bbehsaz/NPDtools-2.3.0-Linux/bin/
#nrpminerpath=/Users/bbehsaz/nrpminer_paper/nrpminer_jan29/nrpminer_april/nrpminer_code/



mkdir -p $outputfolder

python "$nrpminerpath"/scripts/process_correspondence.py "$xaa" $outputfolder

echo $outputfolder"/NRPminer_formatted_correspondence.txt"
if [ -f "$outputfolder"/NRPminer_formatted_correspondence.txt ]; then
	while IFS=$'\t' read -r -a myArray
		do
			genome=$(echo "${myArray}")
			genomefile=$(find "$antismash_folder" -name "*$genome*" -type d -maxdepth 1)
			echo "==============================="
			echo "Processing ... $genomefile"
			python "$nrpminerpath"/nrpminer.py --derepdir "$dereppath" --maxmod 150 --topAA 2 -o "$outputfolder"/NRPDBs/"$genome" --delete_aas "$orfdel" --makeNRPDBonly --antismash_resdir "$genomefile"
			wait 
			#for specname in "${myArray[@]:0}"
			find "$ms_folder"  \( -name "*.mzXML" -o -name "*.mgf" \) -type f -maxdepth 2 | while read specfile;
				do
						specname="$(basename $specfile)"
						echo $specfile
						echo $specname
						find "$outputfolder"/NRPDBs/"$genome" -type d -name "NRPminer_graphs" | while read NRP_DB;
						do 	
							echo $NRPDBs
							python "$nrpminerpath"/nrpminer.py -s "$specfile" --delete_orfs "$orfdel" --derepdir "$dereppath" --maxmod 150 --threads 12 --topAA 2 -o "$outputfolder"/spectrum_nrp_matches/"$specname"_vs_"$genome" --coreNRPDB "$NRP_DB" --coreNRPDBlib "$NRP_DB"/library.info.graphs 
							wait 
						done
						wait
					done
					wait
		done < "$outputfolder"/NRPminer_formatted_correspondence.txt
	wait
fi
