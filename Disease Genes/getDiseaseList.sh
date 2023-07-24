#Obtained browser_source_diseases_summary_CURATED.tsv from DisGeNet website's browse feature, selecting only curated diseases, and then using "Disease class: Mental disorders" and "Semantic type: disease or syndrom" as filters
##The contents of that search wre copied into the file MentalDiseases_DisGeNet.txt
#awk -v FS="\t" -v x="" 'NR>1{x=x"::"$2}END{print x}' browser_source_diseases_summary_CURATED.tsv

#Obtained browser_source_diseases_summary_ANIMAL_MODELS.tsv similarly, but using animal model data
awk -v FS="\t" -v x="" 'NR>1{x=x"::"$2}END{print x}' browser_source_diseases_summary_ANIMAL_MODELS.tsv

#Ids were entered into search and the file C0036341__C0023893__C0005586__C3714756__C0011581__C0001973_disease_gda_summary.tsv was downloaded with the "Download" button.
awk -v FS="\t" -v OFS="\t" 'NR > 1 && $12>.2{gsub(" ", "_", $1); print $3,$1}' C0001849__C0001957__C0002395__C0007384__C0011269__C0014394__C0018614__C0020179__C0021603__C0022336__C0026884__C0027404__C0027821__C0033139__C0035258__C0040517__C0043121__C0155922__C0175691__C0175692__C0206019__C.tsv > diseaseGenes.txt
awk -v FS="\t" -v OFS="\t" 'NR > 1 && $12>.2{gsub(" ", "_", $1); print $3,$1}' C0001306__C0001849__C0002395__C0011263__C0011269__C0020179__C0022336__C0027404__C0035258__C0040517__C0043121__C0175691__C0175692__C0206019__C0242350__C0265338__C0270985__C0276548__C0338451__C0520716__C0795833__C.tsv >> diseaseGenes.txt
sort -k1b,1  diseaseGenes.txt | grep -v -e Variant -e Juvenile -e Onset -e "Familial" -e Sleeper -e Subwakefullness > temp; mv temp diseaseGenes.txt

for file in $(ls /data2/samir/Vmao_Nanopore/GeoMx_Tables/GCvsFLT-SAL/*.mod.txt)
do
	region=$(basename $file | cut -f1 -d"_")
	Region=$(echo $region | awk '{sub("Cortex", "CT"); sub("RCT", "CT"); print}')
	echo $region
	awk -v FS="\t" -v OFS="\t" 'NR>7 && sqrt($4*$4) > .585 && $6 < .05{ print $3,$4}' $file | sort -k1b,1 > ${region}_DEGs.txt
	join -i ${region}_DEGs.txt diseaseGenes.txt | awk -v region=$Region -v OFS="\t" 'BEGIN{print "DEG","brain region","log2 FC", "Associated disease"}{print $1,region,$2,$3}' | uniq | awk -v OFS="\t" '{gsub("_", " "); print}' | awk -v FS="\t" -v OFS="\t" -v gene="" '{if(NR == 1) line=$0; if($1 != gene && NR > 1) {gene=$1; print line; line=$0} else {line=line"; "$4}}' > ${region}_diseaseDEGs.txt
done
cat *_diseaseDEGs.txt | head -1 > Table1.txt
cat *_diseaseDEGs.txt | grep -v DEG >> Table1.txt
cat Table1.txt
echo Total number of disease-associated DEGs:
cut -f1 *_diseaseDEGs.txt | awk '$1 != "DEG"' | sort | uniq | wc -l
