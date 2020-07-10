#!/usr/bin/env python3

import argparse
import yaml
import pprint
from os import makedirs, remove
from sys import stderr
from pathlib import Path
from shutil import copyfile
from subprocess import run
from biopython_helpers import write_bed, make_location, get_gene_data, gene_to_upstream, gene_to_introns 
from metadata import GeneMetaReport


default_gene_ids = [
    #    ACE2        ABO    TMPRSS2
        59272,        28,      7113, # Human
        70008,     80908,     50528, # Mouse
       492331,               494080, # Zebrafish
       712790,    722252,    715138, # Rhesus monkey
    112313373, 112320051, 112306012  # Common vampire bat
    ]


def process_protein_fasta(dataset_dir, dest_dir):
    run(['bash', 'get_proteins.sh', dataset_dir, dest_dir])


def output_location(gene_data, bed_output):
    for gene in gene_data['genes']:
        gene_id = gene['geneId']
        gene_location = make_location(gene['genomicRanges'][0])
        name = f'GeneID_{gene_id}'
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
                name = f'{symbol}_cds:{protein_accession}'
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
                name = f'{symbol}_exon:{transcript_accession}'
                write_bed(bed_output, exon, name)

def output_transcript_variants(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        for transcript in gene['transcripts']:
            transcript_accession = transcript['accessionVersion']
            transcript_name = transcript.get('name', 'transcript').replace(' ', '_')
            transcript_location = make_location(transcript['genomicRange'])
            name = f'{symbol}_{transcript_name}:{transcript_accession}'
            write_bed(bed_output, transcript_location, name)


def output_upstream_regions(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        upstream_collection = gene_to_upstream(gene)
        print('\n●●●●●● BED FORMAT: ●●●●●●')
        for row in upstream_collection:
            name = f'{symbol}_upstream_region:'+','.join(row[1:])
            write_bed(bed_output, row[0], name)


def output_introns(gene_data, bed_output):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        intron_collection = gene_to_introns(gene)
        for row in intron_collection:
            name = f'{symbol}_intron:'+','.join(row[1:])
            write_bed(bed_output, row[0], name)


def process_all(gene_list, output_dir):
    dataset_dir=Path('sars-cov2-gene-data')
    base='sars-cov2-host-genes'
    bed_unsorted = output_dir/(base+'.tmp')
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
    meta_report = GeneMetaReport(gene_list)
    with open(output_dir/(base+'-metadata.tsv'), 'w') as tsv_output:
        meta_report.write_report(tsv_output)

    stderr.write('Done\n')

def main():
    parser = argparse.ArgumentParser(description='Characterization of SARS-CoV-2 host genes.')
    default_gene_ids_string = ', '.join(str(id) for id in default_gene_ids)
    parser.add_argument('--output', nargs='?',
                        default='.',
                        help=f'Output directory. DEFAULT: Current Directory (.)')
    parser.add_argument('genes', nargs='*',
                        default=default_gene_ids,
                        help=f'Input Gene IDs process. DEFAULT: {default_gene_ids_string}')
    args = parser.parse_args()
    output_dir = Path(args.output)
    makedirs(output_dir, exist_ok=True)
    process_all(args.genes, output_dir)


if __name__ == "__main__":
    main()
