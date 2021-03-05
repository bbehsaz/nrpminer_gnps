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
echo $derepdir
echo python "$cyclonovoDir"/scripts/get_by_pepmass.py "$output"_masses.txt "$output"_
cat $allspectfile  | python "$cyclonovoDir"/scripts/get_by_pepmass.py "$output"_masses.txt "$output"_
wait
# cat $allspectfile  | python "$cyclonovoDir"/scripts/get_by_pepmass.py "$massfile" "$output"_

while read line; do if [ ! -d "$output"_"$line"_graphs  ]; then echo "we got there"; "$derepdir"/misc/ad_hoc_python_scripts/graphs_to_lib_info.py "$output"_"$line"_candidateNRPs.GRAPHS.txt "$output"_"$line"_graphs ; fi ; done < "$output"_masses.txt
wait
while read line; do awk '{split($1,a,"/"); print a[length(a)-1]"/"a[length(a)]" "$2" "$3" "$4" "$5}' "$output"_"$line"_graphs/library.info.graphs > temp.txt; mv temp.txt "$output"_"$line"_graphs/library.info.graphs ; done < "$output"_masses.txt
wait


if [ $maxmod != 0 ]; then
	while read line; do 
	echo "BALE?"
	echo $maxmod
	if [ ! -f "$permdir"/"$line"_all_psms.concise.tsv ] || [ ! -d "$permdir"/"$line" ]; then
		# echo "$derepdir"/cycloquest_minimal/varquest.py -o "$outdir"/"$line" --max-charge 1 --max-mod 150 --db-path="$output"_"$line"_graphs "$output"_"$line".mgf -l  "$output"_"$line"_graphs/library.info.graphs; 
		echo "$output"_"$line"_graphs
		"$derepdir"/varquest.py -o "$permdir"/"$line" --max-charge 1 --threads 4 --max-mod="$maxmod" --db-path="$output"_"$line"_graphs --p-value-limit=1e-25 "$output"_"$line".mgf -l  "$output"_"$line"_graphs/library.info.graphs ; 
	
		#One per PSM
		#find "$output"_"$line"_graphs/ -name "com*" | while read gfile; do echo $gfile; echo "$output"_"$line".mgf; "$cyclonovoDir"/scripts/print_score "$output"_"$line".mgf --graph_in $gfile --blind_search --no_merge --no_filter  --concise_output | grep -v "#" | tr " " "_"; done > "$output"_"$line"_all_psms.concise.txt;
		# if [ ++i % 2 == 0 ]; then 
		# 	wait
		# fi
		#wait
		#python "$cyclonovoDir"/scripts/fixconsie.py "$output"_"$line"_all_psms.concise.txt 6 | sort -nk6,6 > "$permdir"/"$line"_all_psms.concise.tsv
		
		wait
	fi
	done < "$output"_masses.txt
else
	echo "Calculating PSMs"
	while read line; do 
	if [ ! -f "$permdir"/"$line"_all_psms.concise.tsv ]; then
		find "$output"_"$line"_graphs/ -name "com*" | while read gfile; do echo $gfile; echo "$output"_"$line".mgf; "$cyclonovoDir"/scripts/print_score "$output"_"$line".mgf --graph_in $gfile  --blind_mod_pos -1 --no_merge --no_filter  --concise_output | grep -v "#" | tr " " "_"; done  > "$output"_"$line"_all_psms.concise.txt;
		# | awk '($0~/^.graph$/){a=$0; getline; b=$0; getline; c=$0; getline; d=$0;getline; e=$0; getline; print a"\t"b"\t"c"\t"d"\t"e"\t"$0}'
		echo "We're done here"
		wait
		python "$cyclonovoDir"/scripts/fixconsie.py "$output"_"$line"_all_psms.concise.txt 3 | sort -nk3,3 | awk '{print $1"\t"$2"\t"$3"\t"$4}' > "$permdir"/"$line"_all_psms.concise.tsv
		wait
	fi;
	# echo "$derepdir"/dereplicator.py -o "$outdir"/"$line" --max-charge 1 --db-path="$output"_"$line"_graphs "$output"_"$line".mgf -l  "$output"_"$line"_graphs/library.info.graphs; 
	# "$derepdir"/dereplicator.py -o "$outdir"/"$line" --max-charge 1 --threads 4 --db-path="$output"_"$line"_graphs  "$output"_"$line".mgf -l  "$output"_"$line"_graphs/library.info.graphs ;
	# python "$cyclonovoDir"/scripts/fixconsie.py "$permdir"/"$line"_all_psms.concise.txt 4 

done < "$output"_masses.txt

fi

wait

rm -R "$outdir"

wait

