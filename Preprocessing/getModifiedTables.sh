for file in $(ls /data2/samir/Vmao_Nanopore/GeoMx_Tables/*/*txt)
do
	outFile=$(echo $file | awk '{sub(".txt", ".mod.txt"); print}')
	awk -v FS="\t" -v OFS="\t" '{if(NR>7 && $NF + 0 != $NF) $NF=2.2; print}' $file > $outFile
done

