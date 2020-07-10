#!/usr/bin/env python3

import argparse
import yaml
import pprint
from os import path,environ, makedirs, remove
from os.path import dirname, realpath
from sys import argv, stderr
from pathlib import Path
from shutil import copyfile
from subprocess import run
from biopython_helpers import write_bed, make_location, get_gene_data, gene_to_upstream, gene_to_introns 
from taxonomy_helpers import get_short_taxname_from_gene
from metadata import GeneMetaReport


default_gene_ids = [
    #    ACE2        ABO    TMPRSS2
        59272,        28,      7113, # Human
        70008,     80908,     50528, # Mouse
       492331,               494080, # Zebrafish
       712790,    722252,    715138, # Rhesus monkey
    112313373, 112320051, 112306012  # Common vampire bat
    ]


def get_script_path():
    return Path(dirname(realpath(argv[0])))


def process_protein_fasta(dataset_dir, dest_dir):
    run(['bash', get_script_path()/'get_proteins.sh', dataset_dir, dest_dir])


def make_bed_name(gene, feature_type=None, rest=None):
    tax_prefix = get_short_taxname_from_gene(gene)
    symbol = gene['symbol']
    gene_id = gene['geneId']
    name = f'{tax_prefix}_{symbol}_GeneID:{gene_id}'
    if feature_type:
        name += '_' + feature_type.replace(' ', '_')
    if rest:
        name += ':' + rest.replace(' ', '_')
    return name


def output_location(gene_data, bed_output):
    for gene in gene_data['genes']:
        gene_location = make_location(gene['genomicRanges'][0])
        name = make_bed_name(gene)
        write_bed(bed_output, gene_location, name)


def output_cds(gene_data, bed_output):
    # TODO: Convert to genomic coordinates.
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        for transcript in gene['transcripts']:
            if not 'cds' in transcript: continue

            cds = transcript['cds']
            protein_accession = transcript['protein']['accessionVersion']
            transcript_location = make_location(transcript['genomicRange'])
            cds_location = make_location(transcript['cds'], transcript_location.strand)
            for cdregion in cds_location.parts:
                name = make_bed_name(gene, 'cds', protein_accession)
                write_bed(bed_output, cdregion, name)


def output_exons(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        for transcript in gene['transcripts']:
            transcript_accession = transcript['accessionVersion']
            # Need strand from genomcRnages, since it's not included in exon locations.
            transcript_location = make_location(transcript['genomicRange'])
            exons_location = make_location(transcript['exons'], transcript_location.strand)
            for exon in exons_location.parts:
                name = make_bed_name(gene, 'exon', transcript_accession)
                write_bed(bed_output, exon, name)

def output_transcript_variants(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        for transcript in gene['transcripts']:
            transcript_accession = transcript['accessionVersion']
            transcript_name = transcript.get('name', 'transcript')
            transcript_location = make_location(transcript['genomicRange'])
            name = make_bed_name(gene, transcript_name, transcript_accession)
            write_bed(bed_output, transcript_location, name)


def output_upstream_regions(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        upstream_collection = gene_to_upstream(gene)

        for row in upstream_collection:
            name = make_bed_name(gene, 'upstream_region', ','.join(row[1:]))
            write_bed(bed_output, row[0], name)


def output_introns(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        intron_collection = gene_to_introns(gene)
        for row in intron_collection:
            name = make_bed_name(gene, 'intron', ','.join(row[1:]))
            write_bed(bed_output, row[0], name)


def process_all(gene_list, output_dir, api_key):
    dataset_dir=Path('sars-cov2-gene-data')
    base='sars-cov2-host-genes'
    bed_unsorted = output_dir/(base+'.tmp')

    selected_genes_text = ', '.join(str(id) for id in gene_list)
    stderr.write(f'Generating SARS-CoV2 host gene report data in: {base}\n')
    stderr.write(f'Selected genes: {selected_genes_text}.\n')
    with open(bed_unsorted, 'w') as bed_output:
        stderr.write('- Fetching gene data\n')
        gene_data = get_gene_data(gene_list, dataset_dir)

        stderr.write('- Processing gene locations\n')
        output_location(gene_data, bed_output)

        stderr.write('- Processing CDS locations\n')
        output_cds(gene_data, bed_output)

        stderr.write('- Processing exons\n')
        output_exons(gene_data, bed_output)

        stderr.write('- Processing transcript variants\n')
        output_transcript_variants(gene_data, bed_output)

        stderr.write('- Processing upstream regions\n')
        output_upstream_regions(gene_data, bed_output)

        stderr.write('- Processing introns\n')
        output_introns(gene_data, bed_output)

    stderr.write('- Processing protien FASTA\n')
    copyfile(dataset_dir/'ncbi_dataset/data/protein.faa', output_dir/(base+'-protein.fasta'))

    stderr.write('- Processing protien domains FASTA\n')
    process_protein_fasta(dataset_dir, output_dir/(base+'-protein-domains.fasta'))

    stderr.write('- Sorting BED output\n')
    run(['sort', '-o', output_dir/(base+'.bed'), '--', bed_unsorted])
    remove(bed_unsorted)

    stderr.write('- Precessing metadata\n')
    meta_report = GeneMetaReport(gene_list, api_key)
    with open(output_dir/(base+'-metadata.tsv'), 'w') as tsv_output:
        meta_report.write_report(tsv_output)

    stderr.write('Done\n')

def main():
    parser = argparse.ArgumentParser(description='Characterization of SARS-CoV-2 host genes.')
    default_gene_ids_string = ', '.join(str(id) for id in default_gene_ids)
    parser.add_argument('--output', nargs='?',
                        default='.',
                        help=f'Output directory. DEFAULT: Current Directory (.)')
    parser.add_argument('--input', nargs='?',
                        help='Input file containing a list of Gene IDs instead of command line')
    parser.add_argument('--api-key', nargs='?',
                        default=environ.get('NCBI_API_KEY', 'ab0568529a7dd0e599fd12b3498f1c8e9e08'),
                        help=f'NCBI API Key for using Entrez faster, with reduced throttling.')
    parser.add_argument('genes', nargs='*',
                        default=default_gene_ids,
                        help=f'Input Gene IDs process. DEFAULT: {default_gene_ids_string}')
    args = parser.parse_args()
    output_dir = Path(args.output)
    makedirs(output_dir, exist_ok=True)
    if args.input:
        if args.genes != default_gene_ids:
            stderr.write('WARNING: Ignoring additional genes listed on command line.\n')    
        gene_set = set()
        with open(args.input) as input:
            for line in input:
                for item in line.split():
                    if item.strip():
                        gene_set.add(int(item))
    else:
        gene_set = set(args.genes)

    if args.api_key:
        environ['NCBI_API_KEY'] = args.api_key

    process_all(gene_set, output_dir, args.api_key)


if __name__ == "__main__":
    main()
