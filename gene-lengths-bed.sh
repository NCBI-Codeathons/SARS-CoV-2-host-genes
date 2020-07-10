#!/bin/sh

##########
#
# Run datasets for list of gene ids that are in a file and process them into a BED file.
# Script looks for geneids.txt
#
##########
     	
for name in $(cat geneids.txt)
do
	../datasets gene-descriptors gene-id $name > scratch.txt
	species=$(cat scratch.txt | jq '.genes[].common_name' -j)
	genename=$(cat scratch.txt | jq '.genes[].common_name, .genes[].description' -j)
	chr=$(cat scratch.txt | jq '.genes[].genomic_ranges[].accession_version' -j)
	begin=$(cat scratch.txt | jq '.genes[].genomic_ranges[].range[].begin' -j)
	begin=$((begin-1)) #convert to 0 counting
	end=$(cat scratch.txt | jq '.genes[].genomic_ranges[].range[].end' -j)
	orientation=$(cat scratch.txt | jq '.genes[].genomic_ranges[].range[].orientation' -j)
	if [ $orientation = "minus" ]; then 
	 	orientation="-"
	 else 
	 	orientation="+"
	 fi
	printf "$chr	$begin	$end	GeneID_$name 1	$orientation\n" >> sortfile.txt
	rm scratch.txt
done
	sort sortfile.txt
	rm sortfile.txt