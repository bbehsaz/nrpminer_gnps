#ms_folder=$1
#antismash_folder_main=$2
#outputfolder=$3
#dereppath=$4
#nrpminerpath=$5
#orfdel=$6
#dereppath=/Users/bbehsaz/NPDtools-2.3.0-Linux/bin/
#nrpminerpath=/Users/bbehsaz/nrpminer_paper/nrpminer_jan29/nrpminer_april/nrpminer_code/



mkdir -p $outputfolder



find $antismash_folder_main -type d -maxdepth 1 | while read antismash_folder;
		do
			echo $antismash_folder
			count=$(ls -1 $antismash_folder/*cluster*.gbk 2>/dev/null | wc -l)
			if [ $count != 0 ]
			then
				genome=$(echo $(basename $antismash_folder))
				echo "==============================="
				echo "Processing ... $genome"
				python "$nrpminerpath"/nrpminer.py --derepdir "$dereppath" --maxmod 150 --topAA 2 -o "$outputfolder"/NRPDBs/"$genome" --delete_aas "$orfdel" --makeNRPDBonly --antismash_resdir "$antismash_folder" > "$outputfolder"/"$genome"_orfDel_"$orfdel".log
				wait 
				
				find "$ms_folder"  \( -name "*.mzXMLt" -o -name "*.mgft" \) -type f -maxdepth 2 | while read specfile;
					do
							specname="$(basename $specfile)"
							# echo $specfile
							echo $specname
							find "$outputfolder"/NRPDBs/"$genome" -type d -name "NRPminer_graphs" | while read NRP_DB;
							do 	
								echo $NRPDBs
								python "$nrpminerpath"/nrpminer.py -s "$specfile" --delete_orfs "$orfdel" --derepdir "$dereppath" --maxmod 150 --threads 12 --topAA 2 -o "$outputfolder"/spectrum_nrp_matches/"$genome"/"$specname" --coreNRPDB "$NRP_DB" --coreNRPDBlib "$NRP_DB"/library.info.graphs 
								wait 
							done
							wait
					done
					wait
			fi
		done
	
