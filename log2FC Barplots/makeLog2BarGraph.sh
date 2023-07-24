echo makeLog2BarGraph.sh
mkdir /data2/samir/Vmao_Nanopore/BarGraph
cd /data2/samir/Vmao_Nanopore/BarGraph
Dir=/home/samir/Vmao_Nanopore


# Pick one Comparison
Compare=GCvsFLT-SAL
echo The comparison is $Compare

for type in CA DG FCortex RCortex
do
	file1=$(ls /data2/samir/Vmao_Nanopore/tables/${type}_GCvsFLT-SAL_top5.txt)
	file=/data2/samir/Vmao_Nanopore/BarGraph/${type}_GCvsFLT-SAL_barplotInput.txt
	sample=$(basename $file _${Compare}_barplotInput.txt)
	awk -v OFS="\t" -v type="CA" 'BEGIN{print type"_gene","Log2","AdjPval"}{print}' $file1 > $file
	Rscript ${Dir}/makeLog2BarGraph.R $file $sample
	echo $type
done
wait
echo done
