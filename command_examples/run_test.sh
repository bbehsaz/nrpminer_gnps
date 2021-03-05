mgf=$1
antismashres=$2
outname=$3
echo $outname
python nrpminer.py -s $mgf --antismash_resdir $antismashres --orfDel 2 --derepdir /Users/bahar/workspace/npd_tools/dereplicator/cycloquest_minimal/bin/ --maxmod 150 --topAA 2 -o test/"$outname" --delete_aas 0 

#for i in {0..2}
#do
# # your-unix-command-here
# echo $i
# python nrpquest_nrpmasses_allderep.py -s $mgf --antismash_resdir $antismashres --delete_orfs $i --derepdir ./../../dereplicator/ --maxmod 150 --topAA 2 -o test/"$outname"/orfs_deleted_"$i" --delete_aas 0
#done

