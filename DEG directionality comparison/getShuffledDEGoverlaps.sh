echo doInterSampleVenns.sh
Dir=/home/samir/Vmao_Nanopore

#Get file of all gene symbols
#cut -f1 /data2/samir/Vmao_Nanopore/heatmaps/CA_GCvsFLT-SAL.allmod > allGenes.txt

#Get up and down DEG files
for region in CA DG FCortex RCortex
do
	GCvsFLTfile=/data2/samir/Vmao_Nanopore/heatmaps/${region}_GCvsFLT-SAL.mod
	SALvsTXfile=/data2/samir/Vmao_Nanopore/heatmaps/${region}_SALvsTX-FLT.mod

	awk '$2 > 0' $GCvsFLTfile > ${region}_GCvsFLTup.txt
	awk '$2 < 0' $GCvsFLTfile > ${region}_GCvsFLTdown.txt
	awk '$2 > 0' $SALvsTXfile > ${region}_SALvsTXup.txt
        awk '$2 < 0' $SALvsTXfile > ${region}_SALvsTXdown.txt
done

#Run the shuffling 1000 times
for i in $(seq 1000)
do

	# This for loop shuffles the file, then makes the up and down direction DEG files from the gene names and their Log2 values made from .mod files made in makeDiffHeatmap.sh.
	for region in CA DG FCortex RCortex
	do
		shuf allGenes.txt > tempshuf &
		SALvsTXupcount=$(cat ${region}_SALvsTXup.txt | wc -l)
                SALvsTXdowncount=$(cat ${region}_SALvsTXdown.txt | wc -l)
		wait

		head -$SALvsTXupcount tempshuf | sort -b > SALvsTXupshuf.txt &
                tail -$SALvsTXdowncount tempshuf | sort -b > SALvsTXdownshuf.txt &
		wait

		join ${region}_GCvsFLTup.txt SALvsTXdownshuf.txt > temp
		join ${region}_GCvsFLTdown.txt SALvsTXupshuf.txt >> temp
		cat temp | wc -l >> ${region}_shuffledOverlaps.txt
	done
done
echo done
