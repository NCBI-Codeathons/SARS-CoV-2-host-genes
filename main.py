#!/usr/bin/env python3

import argparse
import yaml
import pprint
import os
from biopython_helpers import *
from metadata import GeneMetaReport

def process_cds_exons_variants():
    pp = pprint.PrettyPrinter(indent=4)
    # TODO: data_report.yaml is needed for code to work. removed my test file for now
    a_yaml_file = open("data_report.yaml")
    bed_file = open('cds_output.bed', 'w+')
    parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

    parseTranscriptsCds(parsed_yaml_file, pp)
    parseTranscriptsExons(parsed_yaml_file, pp)
    parseTranscriptsName(parsed_yaml_file, pp)
    outputBEDfile(parsed_yaml_file, bed_file, pp)

    a_yaml_file.close()
    bed_file.close()


def parseTranscriptsCds(parsed_yaml_file, pp):
    cds = parsed_yaml_file["genes"][0]["transcripts"][0]["cds"]
    pp.pprint(cds)
    return cds


def parseTranscriptsExons(parsed_yaml_file, pp):
    exons = parsed_yaml_file["genes"][0]["transcripts"][0]["exons"]
    pp.pprint(exons)
    return exons


def parseTranscriptsName(parsed_yaml_file, pp):
    name = parsed_yaml_file["genes"][0]["transcripts"][0]["name"]
    pp.pprint(name)
    return name


def outputBEDfile(parsed_yaml_file, bed_file, pp):
    chromosome = parsed_yaml_file["genes"][0]['chromosomes'][0]
    accessionVersion = parsed_yaml_file["genes"][0]['genomicRanges'][0]['accessionVersion']
    chromStart = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["begin"]
    chromEnd = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["end"]
    orientation = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["orientation"]
    geneId = parsed_yaml_file["genes"][0]["geneId"]
    symbol = parsed_yaml_file["genes"][0]["symbol"]
    name = f"GeneID:{geneId}_{symbol}"
    score = 1

    if orientation == "positive":
        strand = '+'
    elif orientation == "minus":
        strand = '-'
    else:
        strand = '.'

    bed_row = f"{accessionVersion}    {chromStart}    {chromEnd}    {name}    {score}    {strand}\n" 
    pp.pprint(bed_row)
    bed_file.write(f"{accessionVersion}    {chromStart}    {chromEnd}    {name}    {score}    {strand}\n")


def output_upstream_regions(gene_data):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        upstream_collection = gene_to_upstream(gene)
        for row in upstream_collection:
            name = f'{symbol}_upstream_region:'+','.join(row[1:])
            write_bed(row[0], name)

default_gene_ids = [
    #    ACE2        ABO    TMPRSS2
        59272,        28,      7113, # Human
        70008,     80908,     50528, # Mouse
       492331,               494080, # Zebrafish
       712790,    722252,    715138, # Rhesus monkey
    112313373, 112320051, 112306012  # Common vampire bat
    ]

def process_genes(gene_list):
    dest='sars-cov2-gene-data'
    gene_data = get_gene_data(gene_list, dest)
    output_upstream_regions(gene_data)
    os.chdir(f'{dest}/ncbi_dataset/data')
    process_cds_exons_variants()

def main():
    parser = argparse.ArgumentParser(description='Characterization of SARS-CoV-2 host genes.')
    default_gene_ids_string = ', '.join(str(id) for id in default_gene_ids)
    parser.add_argument('genes', nargs='?',
                        default=default_gene_ids,
                        help=f'Input Gene IDs process. DEFAULT: {default_gene_ids_string}')
    args = parser.parse_args()
    install_datasets()
    process_genes(args.genes)
    meta_report = GeneMetaReport(args.genes)
    meta_report.write_report()


if __name__ == "__main__":
    main()
