allreconstfile=$1

allspectfile=$2

output=$3

cyclonovoDir=$4

threshcore=2
err=$5

cat $allreconstfile | awk '($(NF)>7){print $(NF-3)}' | sort | uniq > "$output"_masses.txt

echo "Calculating P-values and cleaning up ... "
echo "BALE?"
while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>5){n++; if (n<10000){print $1}}' | sed 's/,/ /g' > "$output"_"$line"_reconst.txt  ; done < "$output"_masses.txt
while read line; do grep "$line" "$allreconstfile" | awk '(NF>3 && $(NF-2)>5){n++; if (n<10000){print $0}}' > "$output"_"$line"_allinfo.txt  ; done < "$output"_masses.txt

echo "$output"_masses.txt


cat $allspectfile  | python "$cyclonovoDir"/scripts/get_by_pepmass.py "$output"_masses.txt "$output"_


while read line; do if [ -f "$output"_"$line".mgf ]; then "$cyclonovoDir"/scripts/print_score "$output"_"$line".mgf "$output"_"$line"_reconst.txt --mass_seq_in --no_filter --no_merge --make_cyclic --multiple_seq_file --concise_output | grep -v "^#" > "$output"_"$line"_pvals.txt; fi; done < "$output"_masses.txt 

while read line; do if [ -f "$output"_"$line".mgf ]; then paste "$output"_"$line"_allinfo.txt "$output"_"$line"_pvals.txt | awk '{print $(NF-1),$0}' | sort -nr | cut -f2- -d' ' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  > "$output"_"$line"_reconst_pvals.txt; fi; done < "$output"_masses.txt 

while read line; do if [ -f "$output"_"$line"_allinfo.txt ]; then rm "$output"_"$line"_allinfo.txt; fi; done < "$output"_masses.txt

while read line; do if [ -f "$output"_"$line"_pvals.txt ]; then rm "$output"_"$line"_pvals.txt; fi; done < "$output"_masses.txt

while read line; do if [ -f  "$output"_"$line".mgf ]; then rm  "$output"_"$line".mgf; fi; done < "$output"_masses.txt

while read line; do if [ -f  "$output"_"$line"_reconst.txt ]; then rm  "$output"_"$line"_reconst.txt; fi; done < "$output"_masses.txt

rm $allreconstfile
while read line; do if [ -f  "$output"_"$line"_reconst_pvals.txt ]; then cat "$output"_"$line"_reconst_pvals.txt >> "$allreconstfile".final.txt ;fi; done < "$output"_masses.txt
while read line; do if [ -f  "$output"_"$line"_reconst_pvals.txt ]; then head -n1 "$output"_"$line"_reconst_pvals.txt >> summary.txt ;fi; done < "$output"_masses.txt
while read line; do if [ -f  "$output"_"$line"_reconst_pvals.txt ]; then rm  "$output"_"$line"_reconst_pvals.txt; fi; done < "$output"_masses.txt
