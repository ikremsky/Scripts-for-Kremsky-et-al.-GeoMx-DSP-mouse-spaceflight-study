echo countDEgenes.sh
cd /data2/samir/Vmao_Nanopore/tables
Dir=/home/samir/Vmao_Nanopore

# Pick one Comparison
Compare1=GCvsFLT-SAL
Compare2=SALvsTX-FLT
Compare3=SALvsTX-GC
echo The comparison is $Compare1 and $Compare2

# Makes a table of the significant genes (q < .05 and |log FC| > .585, equivalant to FC magnitude of atleast 150%)
for Comp in $Compare1 $Compare2 $Compare3
do
	for file in $(ls ../GeoMx_Tables/${Comp}/*_${Comp}.mod.txt)
	do
	        sample=$(basename $file _${Comp}.mod.txt)
	        awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{gsub("All Targets","All_Targets");print $3}' $file | sort -b | uniq > ${sample}_${Comp}_DEgenes.txt
	done
done
wait
wc -l *GCvsFLT-SAL*_DEgenes.txt

for sample in CA DG FCortex RCortex
do
	echo $sample
	join ${sample}_${Compare1}_DEgenes.txt ${sample}_${Compare2}_DEgenes.txt > ${sample}_joinedDEgenes.txt
done
wait
wc -l *joinedDEgenes.txt
head *joinedDEgenes.txt
echo done
