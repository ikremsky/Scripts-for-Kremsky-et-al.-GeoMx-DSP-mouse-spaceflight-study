echo makeDiffHeatmap.sh
mkdir /data2/samir/Vmao_Nanopore/heatmaps
cd /data2/samir/Vmao_Nanopore/heatmaps
Dir=/home/samir/Vmao_Nanopore

# Pick one comparison
Compare1=GCvsFLT-SAL
Compare2=SALvsTX-FLT
Compare3=SALvsTX-GC
Compare4==GCvsFLT-TX
Compare=SALvsTX
Compare=${Compare1}_${Compare2}

rm *join* genes_*.txt

# This script makes diffheatmaps for FCortex, Cortex, CA, and DG samples FOR GCvsFLT-SAL SALvsTX-FLT. Must run makeDiffHeatmap.sh first to get .mod files

for file in $(ls ../GeoMx_Tables/*/*.mod.txt)
do
        awk -v FS="\t" -v OFS="\t" 'NR>7{gsub("All Targets","All_Targets");print $3,$4}' $file | sort -k1b,1 > $(basename $file .mod.txt).allmod &
done

cut -f1 *${Compare1}.mod | sort -b | uniq > genes_${Compare1}.txt
wait
head genes_${Compare1}.txt

for Comp in $Compare1  $Compare2 $Compare3 $Compare4
do
	for file in $(ls *_${Comp}.allmod)
	do
		name=$(basename $file .allmod)
		join genes_${Compare1}.txt $file | awk -v name=$name -v OFS="\t" 'BEGIN{print "gene",name}{print $1,$2}' > ${name}.joined &
	done
done
wait

##paste *${Compare}*.joined | awk -v Compare=$Compare -v OFS="\t" '{gsub("_"Compare,""); gsub("RCortex","Cortex");print $1,$2,$4,$6,$8}' > final.table_${Compare}.txt
paste *.joined | awk -v OFS="\t" '{gsub("GCvsFLT","GCvFLT"); gsub("SALvsTX","SALvBuOE"); gsub("Cortex","CT"); gsub("RCT", "CT"); print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24}' > final.table_${Compare}.txt
Rscript ${Dir}/makeDiffHeatmap2.R final.table_${Compare}.txt $Compare
echo done
