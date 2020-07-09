#!/bin/bash

datasets_dir=$1
output_dir=$2

usage() {
        echo "Usage: $0 <DATASETS-DIR> <OUTPUT-DIR>"
        exit 1
}

[[ $# == 2 ]] || usage

protein_fasta=$datasets_dir/ncbi_dataset/data/protein.faa
[[ -r $protein_fasta ]] || {
       echo >&2 "ERROR: Must provide input to NCBI Datasets gene data with protein.faa"
       exit 1
}

[[ -d $output_dir ]] || {
       echo >&2 "ERROR: Must provide output directory"
       exit 1
}

echo -n "" > tempRegion

for i in $(grep '>' "$datasets_dir/ncbi_dataset/data/protein.faa" | cut -f 1 -d ' '| sed 's/>//')
do
        echo $i
        esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd sig_peptide   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "sig_peptide\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd mat_peptide   product INSDInterval_from INSDInterval_to sub_sequence | awk '{print "mat_peptide\t"$0}' >> tempRegion
done

paste tempRegion <(cat tempRegion| cut -f 3 | cut -f 1 -d '/') | awk  'BEGIN {FS="\t"}{print ">"$2"|"$1"|"$4"-"$5"|"$8"\n"$6}' > "$output_dir/feature.fasta"
