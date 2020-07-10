#!/bin/bash

set -o errexit
set -o pipefail

datasets_dir=$1
output=$2

usage() {
        echo "Usage: $0 <DATASETS-DIR> <OUTPUT>"
        exit 1
}

[[ $# == 2 ]] || usage

protein_fasta=$datasets_dir/ncbi_dataset/data/protein.faa
[[ -r $protein_fasta ]] || {
       echo >&2 "ERROR: Must provide input to NCBI Datasets gene data with protein.faa"
       exit 1
}

touch -- "$output" || {
       echo >&2 "ERROR: Could not create output: $output"
       exit 1
}

echo -n "" > tempRegion

for i in $(grep '>' "$datasets_dir/ncbi_dataset/data/protein.faa" | cut -f 1 -d ' '| sed 's/>//')
do
        echo >&2 "  - Processing FASTA for protein: $i"
        esearch -db protein -query $i |efetch -format gpc|xtract -insd Region   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "Region\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd sig_peptide   region_name INSDInterval_from INSDInterval_to sub_sequence | awk '{print "sig_peptide\t"$0}' >> tempRegion
        esearch -db protein -query $i |efetch -format gpc|xtract -insd mat_peptide   product INSDInterval_from INSDInterval_to sub_sequence | awk '{print "mat_peptide\t"$0}' >> tempRegion
done

paste tempRegion <(cat tempRegion| cut -f 3 | cut -f 1 -d '/') | awk  'BEGIN {FS="\t"}{print ">"$2"|"$1"|"$4"-"$5"|"$8"\n"$6}' > "$output"
