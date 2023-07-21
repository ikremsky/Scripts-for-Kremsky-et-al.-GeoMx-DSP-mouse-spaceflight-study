echo makeDiffHeatmap.sh
mkdir /data2/samir/Vmao_Nanopore/heatmaps
cd /data2/samir/Vmao_Nanopore/heatmaps
Dir=/home/samir/Vmao_Nanopore

#rm *mod*
# Pick one Comparearison
#Compare=GCvsFLT-SAL
#Compare=SALvsTX-FLT
#Compare=SALvsTX-GC
Compare=GCvsFLT-TX

#rm *mod* *join* *final* genes_*.txt

# This script makes diffheatmaps for FCortex, Cortex, CA, and DG samples
for file in $(ls ../GeoMx_Tables/*${Compare}/*_${Compare}.mod.txt)
do
	awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{gsub("All Targets","All_Targets");print $3,$4}' $file | sort -k1b,1 > $(basename $file .mod.txt).mod &
done
wait

for file in $(ls ../GeoMx_Tables/*${Compare}/*_${Compare}.mod.txt)
do
        awk -v FS="\t" -v OFS="\t" 'NR>7{gsub("All Targets","All_Targets");print $3,$4}' $file | sort -k1b,1 > $(basename $file .mod.txt).allmod &
done
wait

cut -f1 *${Compare}.mod | sort -b | uniq > genes_${Compare}.txt




for file in $(ls *_${Compare}.allmod)
do
	name=$(basename $file .allmod)
	join genes_${Compare}.txt $file | awk -v name=$name -v OFS="\t" 'BEGIN{print "gene",name}{print $1,$2}' > ${name}.joined &
done
wait




#for file in $(ls *_${Compare}.mod)
#do
#	name=$(basename $file .mod)
#	join -v 1 genes_${Compare}.txt $file | awk -v OFS="\t" '{print $1,"0"}' > $(basename $file .mod).joined
#	cat $file $(basename $file .mod).joined | sort -k1b,1 | awk -v name=$name -v OFS="\t" 'BEGIN{print "Gene", name}{print $1,$2}' > $(basename $file .mod).final.joined
#done
#wait

paste *${Compare}.joined | awk -v Compare=$Compare -v OFS="\t" '{gsub("_"Compare,""); print $1,$2,$4,$6,$8}' > final.table_${Compare}.txt

Rscript ${Dir}/makeDiffHeatmap.R final.table_${Compare}.txt $Compare
echo done
