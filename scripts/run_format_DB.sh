#this script writes the putative NRP database onto the disk, runs VarQuest/Dereplicator to calcualte the PSMS. 
output=$1
allspectfile=$2
cyclonovoDir=$3
threshcore=4
outdir=$4
derepdir=$5
err=$6
tempdir=temp
maxmod=$7
permdir=$8
groupname=$9
threads=${10}

echo "Core NRP database generated ..." "$output"_graphs
ls -lth "$output"_candidateNRPs.GRAPHS.txt
echo "Writing the core NRP database onto the disk ..." 
python "$cyclonovoDir"/scripts/graphs_to_lib_info_nrpminer_version.py "$output"_candidateNRPs.GRAPHS.txt "$output"_graphs
awk '{split($1,a,"/"); print a[length(a)-1]"/"a[length(a)]" "$2" "$3" "$4" "$5}' "$output"_graphs/library.info.graphs > temp.txt; mv temp.txt "$output"_graphs/library.info.graphs 
wait

# echo "Calculating coreNRP-Spectrum-Matches ... "
# if [ $maxmod != 0 ]; then
# 	# "$derepdir"/varquest.py -o "$permdir"/psms_"$groupname" --max-charge 1 --threads "$threads" --db-path="$output"_graphs --max-mod="$maxmod" --p-value-limit=1e-15 "$allspectfile" -l  "$output"_graphs/library.info.graphs
# 	# "$derepdir"/varquest.py -o "$permdir"/psms_"$groupname" --max-charge 1 --threads "$threads" --db-path="$output"_graphs --max-mod="$maxmod" --p-value-limit=1e-15 "$allspectfile" -l  "$output"_graphs/library.info.graphs --fdr-limit=5
# 	"$derepdir"/varquest.py -o "$permdir"/psms_"$groupname" --max-charge 1 --threads "$threads" --db-path="$output"_graphs --max-mod="$maxmod" --p-value-limit=1e-15 "$allspectfile" -l  "$output"_graphs/library.info.graphs --fdr
# else
# 	"$derepdir"/dereplicator.py -o "$permdir"/psms_"$groupname"_no_modifications --max-charge 2 --threads "$threads" --db-path="$output"_graphs --max-mod="$maxmod" --p-value-limit=1e-15 "$allspectfile" -l  "$output"_graphs/library.info.graphs --fdr
# fi
# wait

# echo "Saving the NRPminer core NPR database to "$outdir
# zip -qr "$outdir"/zipped_nrpminer_database.zip "$output"_graphs
# # rm -R "$output"_graphs
# wait

