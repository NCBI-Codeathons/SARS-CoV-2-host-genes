#!/bin/sh

##########
#
# Run datasets for list of gene ids that are in a file.
#
##########
     	
for name in $(cat geneids.txt)
do
	species=$(../datasets gene-descriptors gene-id $name | jq '.genes[].common_name' -j)
	genename=$(../datasets gene-descriptors gene-id $name | jq '.genes[].common_name, .genes[].description' -j)
	chr=$(../datasets gene-descriptors gene-id $name | jq '.genes[].genomic_ranges[].accession_version' -j)
	begin=$(../datasets gene-descriptors gene-id $name | jq '.genes[].genomic_ranges[].range[].begin' -j)
	begin=$((begin-1)) #convert to 0 counting
	end=$(../datasets gene-descriptors gene-id $name | jq '.genes[].genomic_ranges[].range[].end' -j)
	orientation=$(../datasets gene-descriptors gene-id $name | jq '.genes[].genomic_ranges[].range[].orientation' -j)
	if [ $orientation = "minus" ]; then 
	 	orientation="-"
	 else 
	 	orientation="+"
	 fi
	printf "$chr	$begin	$end	GeneID_$name 1	$orientation\n"
	done