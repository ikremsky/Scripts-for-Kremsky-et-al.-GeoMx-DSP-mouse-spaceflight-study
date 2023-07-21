awk -v OFS="\t" '{gsub("Cortex", "CT"); gsub("RCT", "CT"); print}' /data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt > GeneExpressionTable.txt
Rscript PCA_allSamples.R
