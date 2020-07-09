
import yaml
from biopython_helpers import *

def main():
    install_datasets()
    gene_data = get_gene_data([59272])

    output_upstream_regions(gene_data)

    parseTranscriptsCds(gene_data)
    parseTranscriptsExons(gene_data)
    parseTranscriptsName(gene_data)


def parseTranscriptsCds(gene_data):
    cds = gene_data["genes"][0]["transcripts"][0]["cds"]
    print(f'●●●●●● CDS: ●●●●●●\n{cds}\n')
    return cds


def parseTranscriptsExons(gene_data):
    exons = gene_data["genes"][0]["transcripts"][0]["exons"]
    print(f'●●●●●● EXONS: ●●●●●●\n{exons}\n')
    return exons


def parseTranscriptsName(gene_data):
    name = gene_data["genes"][0]["transcripts"][0]["name"]
    print(f'●●●●●● NAME: ●●●●●●\n{name}\n')
    return name


def output_upstream_regions(gene_data):
    for gene in gene_data['genes']:
        symbol = gene["symbol"]
        upstream_collection = gene_to_upstream(gene)
        print('\n●●●●●● BED FORMAT: ●●●●●●')
        for row in upstream_collection:
            name = f'{symbol}_upstream_region:'+','.join(row[1:])
            write_bed(row[0], name)


if __name__ == "__main__":
    main()
