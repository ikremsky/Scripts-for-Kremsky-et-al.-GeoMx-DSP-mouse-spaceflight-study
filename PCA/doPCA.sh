echo doPCA.sh
mkdir /data2/samir/Vmao_Nanopore/PCA
cd /data2/samir/Vmao_Nanopore/PCA

x=""
for celltype in CA DG FCortex RCortex
do
	x=$(echo $x" "${celltype}_GC_SAL ${celltype}_FLT_SAL)
#	echo $x
#	Rscript /home/samir/Vmao_Nanopore/PCA.R ${celltype} ${celltype}_GC_SAL ${celltype}_FLT_SAL &
done
wait
echo $x
Rscript /home/samir/Vmao_Nanopore/PCA_All_Reps.R $x
echo done
