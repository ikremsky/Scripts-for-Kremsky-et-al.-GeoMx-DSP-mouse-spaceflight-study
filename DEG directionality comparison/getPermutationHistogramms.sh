echo getShuffledBarplots.sh

#Files for runnign this script were generated in BarPlotsFromVennData.sh and getShuffledDEGoverlaps.sh

#Running the R script using the format: Title, Overlaps 1-4 (DD,DU,UD,UU), Total DE gene outputs 1-4 (GCvsFLT down then up, SALvsTX down then up)

Rscript histogram_shuffled.R CA $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/CA_GCvsFLT_*_vs_*) $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/CA_*.*) CA_shuffledOverlaps.txt

Rscript histogram_shuffled.R DG $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/DG_GCvsFLT_*_vs_*) $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/DG_*.*) DG_shuffledOverlaps.txt

Rscript histogram_shuffled.R FCT $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/FCortex_GCvsFLT_*_vs_*) $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/FCortex_*.*) FCortex_shuffledOverlaps.txt

Rscript histogram_shuffled.R CT $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/RCortex_GCvsFLT_*_vs_*) $(ls /home/sjustinen/V_Mao/VennDiagram/Counts/RCortex_*.*) RCortex_shuffledOverlaps.txt

echo done
