echo makeNormalizedReadCountHeatmap.sh
mkdir /data2/samir/Vmao_Nanopore/heatmaps
cd /data2/samir/Vmao_Nanopore/heatmaps
Dir=/home/samir/Vmao_Nanopore

# The commands below pull the rows of the Log2 values that match with the Normalized Readcount Gene Expression table
# Hclust_Log2Table.txt file made in /home/samir/Vmao_Nanopore/makeDiffHeatmap.sh script
# GeneExpressionTable.txt file made in /home/samir/Vmao_Nanopore/makeGeneExpressionTable.sh script

awk 'NR>1{print$1}' Hclust_Log2Table.txt | sort -k1b,1 > Genes
sort -k1b,1 ../tables/GeneExpressionTable.txt > GeneExpressionTable.sorted.txt
head -1 ../tables/GeneExpressionTable.txt > GeneExpressionTable.sorted.mod
join Genes GeneExpressionTable.sorted.txt | awk -v OFS="\t" '{gsub(" ","\t");print}' >> GeneExpressionTable.sorted.mod
Rscript ${Dir}/makeNormalizedReadCountHeatmap.R GeneExpressionTable.sorted.mod Hclust_Log2Table.txt

x=""
for celltype in CA DG FCortex RCortex
do
        x=$(echo $x" "${celltype})
done
wait
Rscript ${Dir}/makeNormalizedReadCountHeatmap_automated.R $x
echo done
