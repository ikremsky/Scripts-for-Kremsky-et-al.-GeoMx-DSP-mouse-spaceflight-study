echo makeGeneExpressionTable.sh
mkdir /data2/samir/Vmao_Nanopore/tables
cd /data2/samir/Vmao_Nanopore/tables
#This script made the input data to PCA.R for all samples
# ${sample}Table.txt obtained from Geo DSP website, from CA,DG,FCT,RCT excel sheet TargetCountMatrix
for sample in CA RCortex DG FCortex
do
	awk -v sample=$sample -v OFS="\t" '{gsub("Full ROI",sample);print}' ${sample}Table.txt | awk -v OFS="\t" '{gsub(/ \| /,"_");print}' > ${sample}Table.temp.txt
	head -1 ${sample}Table.temp.txt | awk '{gsub("\t","\n");print}' | awk -v FS="_" -v OFS="\t" '{print $5"_"$1"_"$2"_"$4}' | paste -d "\t" -s | awk -v OFS="\t" '{gsub("_TargetName__","TargetName");print}' > ${sample}Table.mod.txt
	awk -v OFS="\t" 'NR>1{print}' ${sample}Table.temp.txt | sort >> ${sample}Table.mod.txt
	cut -f 2- ${sample}Table.mod.txt >  ${sample}.Values
done
wait
cut -f1 CATable.mod.txt > Genes
paste Genes *.Values > GeneExpressionTable.txt

for sample in CA RCortex DG FCortex
do
        paste Genes ${sample}*.Values > ${sample}_GeneExpressionTable.txt
done
wait
#rm *.Values *.mod.txt *.temp.txt
echo done
