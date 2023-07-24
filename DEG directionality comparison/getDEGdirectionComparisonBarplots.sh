#Stephen Justinen
echo start
cd ./VennDiagram/Counts

#Generating files containing only the overlapping DE genes, comparing the GCvsFLT to SALvsBuOE Down/Up results

for file in $(ls *_File_GCvsFLT*)
do
	x=1
	for file2 in $(ls *$(echo $file | cut -d'_' -f 1)*_File_SALvsTX* )
	do
		comm -12 $file $file2 > $file$x
		x=$(($x+1))
	done
done

#Renaming Files generated above, removing the numberic extension and replacing it with the full comparision for clairity when reviewing files.

for file in $(ls *Down1)
do
	 mv $file $(echo $file | cut -d'_' -f 1)_GCvsFLT_Down_vs_SALvsBuOE_Down
done

for file in $(ls *Down2)
do
         mv $file $(echo $file | cut -d'_' -f 1)_GCvsFLT_Down_vs_SALvsBuOE_UP
done

for file in $(ls *Up1)
do
         mv $file $(echo $file | cut -d'_' -f 1)_GCvsFLT_UP_vs_SALvsBuOE_Down
done

for file in $(ls *Up2)
do
         mv $file $(echo $file | cut -d'_' -f 1)_GCvsFLT_UP_vs_SALvsBuOE_UP
done
wait
cd Graphs

#Running the R script using the format: Title, Overlaps 1-4 (DD,DU,UD,UU), Total DE gene outputs 1-4 (GCvsFLT down then up, SALvsTX down then up)

Rscript ../DEGdirectionComparisonBarplots.R CA $(ls ../CA_GCvsFLT_*_vs_*) $(ls ../CA_*.*)

Rscript ../DEGdirectionComparisonBarplots.R DG $(ls ../DG_GCvsFLT_*_vs_*) $(ls ../DG_*.*)

Rscript ../DEGdirectionComparisonBarplots.R FCT $(ls ../FCortex_GCvsFLT_*_vs_*) $(ls ../FCortex_*.*)

Rscript ../DEGdirectionComparisonBarplots.R CT $(ls ../RCortex_GCvsFLT_*_vs_*) $(ls ../RCortex_*.*)

echo done
