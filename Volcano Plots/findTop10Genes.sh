echo findTop5Genes.sh
cd /data2/samir/Vmao_Nanopore/tables
Dir=/home/samir/Vmao_Nanopore

# Pick one Comparison
Compare=GCvsFLT-SAL
#Compare=SALvsTX-FLT
#Compare=SALvsTX-GC
echo The comparison is $Compare

# Makes a table with the gene, Log2, and -log10 adjusted pvalue with -log10 adjusted pvalue > 1.3 for each sample
for file in $(ls ../GeoMx_Tables/${Compare}/*_${Compare}.mod.txt)
do
	sample=$(basename $file _${Compare}.txt)
	awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{gsub("All Targets","All_Targets");print $3,$4,$NF,$4}' $file | awk -v OFS="\t" '{gsub("-","",$2);print}' | sort -k2nr,2 | head |  awk -v FS="\t" -v OFS="\t" '{print $1,$4,$3}' | sort -k1b,1 > $(basename $file .mod.txt)_top10.txt &
	awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{gsub("All Targets","All_Targets");print $3,$4,$NF,$4}' $file | awk -v OFS="\t" '{gsub("-","",$2);print}' | sort -k2nr,2 | head -5 |  awk -v FS="\t" -v OFS="\t" '{print $1,$4,$3}' > $(basename $file .mod.txt)_top5.txt &
	echo $sample
	awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{gsub("All Targets","All_Targets");print $3,$4,$NF,$4}' $file | awk -v OFS="\t" '{gsub("-","",$2);print}' | sort -k2nr,2 | head |  awk -v FS="\t" -v OFS="\t" '{print $1,$4,$3}' | head -5 | awk -v x="" '{x=x", "$1}END{sub(", ", "", x); print x}'
done
wait
