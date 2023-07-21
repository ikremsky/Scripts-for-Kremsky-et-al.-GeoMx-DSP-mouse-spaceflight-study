echo getVolcanoPlot.sh
echo makeDiffHeatmap.sh
mkdir /data2/samir/Vmao_Nanopore/VolcanoPlots
cd /data2/samir/Vmao_Nanopore/VolcanoPlots
Dir=/home/samir/Vmao_Nanopore

# Pick one Comparearison
Compare=GCvsFLT-SAL


# The comand below is for making volcano plots from the GeoMx DSP data.
for file in $(ls ../GeoMx_Tables/*${Compare}/*_${Compare}.mod.txt)
do
	sample=$(basename $file _${Compare}.mod.txt)
	echo $sample
	awk -v FS="\t" -v OFS="\t" 'NR>6{gsub("Adjusted pvalue","Adjusted_pvalue"); gsub("Target name","gene_symbol",$3); gsub("Log2","log2FoldChange") ;gsub("-log10 adjusted pvalue","NeglogPadjValue",$NF); print $3,$4,$6,$NF}' $file > ${sample}_input_VolcanoPlot.txt

	awk -v OFS="\t" '{print}' ../tables/${sample}_${Compare}_top10.txt > ${sample}_top10Genes_VolcanoPlot.txt

	Rscript /home/samir/Vmao_Nanopore/getVolcanoPlot_ggplot.R ${sample}_input_VolcanoPlot.txt $sample ${sample}_top10Genes_VolcanoPlot.txt
done
wait
echo done
