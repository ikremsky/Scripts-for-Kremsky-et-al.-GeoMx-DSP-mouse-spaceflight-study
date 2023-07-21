#Pathways affected by oxidative stress: VEGF signaling (CA), cAMP signaling pathway (DG), FoxO signaling pathway (FCT), ECM-receptor interaction (CT)

#python getKEGGpathwayID.py > pathways.txt

Compare1=GCvsFLT-SAL
Compare2=SALvsTX-FLT

for pathway in "ECM-receptor interaction RCortex" "FoxO signaling pathway FCortex" "Calcium signaling pathway DG" "Pyruvate metabolism CA"
do
	name=$(echo $pathway| cut -f1 -d" ")
	region=$(echo $pathway| awk '{print $NF}')
	echo $region
	Pathway=$(echo $pathway| awk '{$NF=""; print}')
	id=$(grep "$Pathway" pathways.txt | cut -f1)
	echo $Pathway $id
#	python getKEGGgenes.py $id | sort -b > ${name}_genes.txt
	wc -l ${name}_genes.txt
	join ${name}_genes.txt /data2/samir/Vmao_Nanopore/tables/${region}_GCvsFLT-SAL_DEgenes.txt > ${region}_${name}.txt
	join ${name}_genes.txt /data2/samir/Vmao_Nanopore/tables/${region}_SALvsTX-FLT_DEgenes.txt >> ${region}_${name}.txt
	sort -b ${region}_${name}.txt | uniq > temp; mv temp ${region}_${name}.txt
	wc -l ${region}_${name}.txt
	for Comp in $Compare1  $Compare2
	do
		CompF=$(echo $Comp | awk '{sub("TX", "BuOE"); print}')
		file=/data2/samir/Vmao_Nanopore/heatmaps/${region}_${Comp}.allmod
#               join ${region}_${name}.txt $file | awk -v name=$CompF -v OFS="\t" 'BEGIN{print "gene",name}{print $1,$2}' > ${region}_${name}_${Comp}.joined &
		join ${name}_genes.txt $file | awk -v name=$CompF -v OFS="\t" 'BEGIN{print "gene",name}{print $1,$2}' > ${region}_${name}_${Comp}.joined &
	done
	wait

	paste ${region}_${name}*.joined | head -1 | awk -v OFS="\t" '{gsub("GCvsFLT","GCvFLT"); gsub("SALvsBuOE","SALvBuOE"); print $1,$2,$4}' > final.table_${region}_${name}.txt
	paste ${region}_${name}*.joined | awk 'NR > 1' | awk -v OFS="\t" '{if(sqrt($2*$2) > sqrt($4*$4)) x=sqrt($2*$2); else x=sqrt($4*$4); print $1,$2,$4,x}' | sort -k4nr,4 | cut -f1-3 | head -20 >> final.table_${region}_${name}.txt
	Rscript makeDiffHeatmap2.R final.table_${region}_${name}.txt ${region}_${name}
done


#Heatmap of Alzheimer's disease genes
grep Alzheimer diseaseGenes.txt | cut -f1 | sort -b | uniq >  diseaseGenes_Alzheimer.txt
for region in CA DG FCortex RCortex
do
	name=Alzheimer
        for Comp in $Compare1  $Compare2
        do
		CompF=$(echo $Comp | awk '{sub("TX", "BuOE"); print}')
                file=/data2/samir/Vmao_Nanopore/heatmaps/${region}_${Comp}.allmod
                join -i $file diseaseGenes_Alzheimer.txt | awk -v name=$CompF -v OFS="\t" 'BEGIN{print "gene",name}{print $1,$2}' > ${region}_${name}_${Comp}.joined &
        done
        wait
	paste ${region}_${name}*.joined | head -1 | awk -v OFS="\t" '{gsub("GCvsFLT","GCvFLT"); gsub("SALvsBuOE","SALvBuOE"); print $1,$2,$4}' > final.table_${region}_${name}.txt
	paste ${region}_${name}*.joined | awk 'NR > 1' | awk -v OFS="\t" '{if(sqrt($2*$2) > sqrt($4*$4)) x=sqrt($2*$2); else x=sqrt($4*$4); print $1,$2,$4,x}' | sort -k4nr,4 | cut -f1-3 | head -20 >> final.table_${region}_${name}.txt
	echo $region
	Rscript makeDiffHeatmap2.R final.table_${region}_${name}.txt ${region}_${name}
done
