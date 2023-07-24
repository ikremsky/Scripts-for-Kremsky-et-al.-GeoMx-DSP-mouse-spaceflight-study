#Stephen Justinen
echo start
#/data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt

#cut -f2-13 -d$'\t' /data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt > brain_1.txt
#cut -f14-25 -d$'\t' /data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt > brain_2.txt
#cut -f26-37 -d$'\t' /data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt > brain_3.txt
#cut -f38-49 -d$'\t' /data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt > brain_4.txt

Rscript correlateReplicates_brain.R brain_1.txt CA
Rscript correlateReplicates_brain.R brain_2.txt DG
Rscript correlateReplicates_brain.R brain_3.txt FCT
Rscript correlateReplicates_brain.R brain_4.txt CT

#Reminder for R Code:
#paste temp1 temp2 | awk '$1 > 1 || $2 > 1' > R_Graphs/${type}_${file}.expressed

echo done

montage -tile 2x2 -geometry 600x600 -density 600 CA*.png Correlation_Heatmaps_CA_montage.png
montage -tile 2x2 -geometry 600x600 -density 600 DG*.png Correlation_Heatmaps_DG_montage.png
montage -tile 2x2 -geometry 600x600 -density 600 FCT*.png Correlation_Heatmaps_FCT_montage.png
montage -tile 2x2 -geometry 600x600 -density 600 CT*.png Correlation_Heatmaps_CT_montage.png
