echo doInterSampleVenns.sh
cd ./VennDiagram
Dir=/home/samir/Vmao_Nanopore

# This for loop makes the up and down direction files from the gene names and their Log2 values made from .mod files made in makeDiffHeatmap.sh.
for file in $(ls /data2/samir/Vmao_Nanopore/heatmaps/*mod*)
do
        name=$(basename $file .mod | awk '{sub("-SAL", ""); sub("-FLT", ""); print}')
        awk -v OFS="\t" '{if($2>0) print}' $file > ${name}.UpReg
        awk -v OFS="\t" '{if($2<0) print}' $file > ${name}.DownReg
done
echo done
