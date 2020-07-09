#!/usr/bin/env python3

import sys
import zipfile
import yaml
import requests
import ncbi.datasets
from os import stat, chmod, remove
from shutil import copyfileobj, rmtree
from stat import S_IEXEC
from subprocess import run
from Bio.SeqFeature import FeatureLocation, CompoundLocation



def install_datasets():
    path = 'datasets'
    url = 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
    r = requests.get(url)
    open(path, 'wb').write(r.content)
    st = stat(path)
    chmod(path, st.st_mode | S_IEXEC)


def download_gene_data_via_command_line(gene_list, dest):
    run(['./datasets', 'download', 'gene', ','.join(str(gene) for gene in gene_list)])
    run(['unzip', '-qod', dest, 'ncbi_dataset.zip'])
    remove('ncbi_dataset.zip')


def download_gene_data_via_api(gene_ids, dest):
    with ncbi.datasets.ApiClient() as api_client:
        api_instance = ncbi.datasets.DownloadApi(api_client)
        include_sequence_type = ['SEQ_TYPE_GENE', 'SEQ_TYPE_RNA', 'SEQ_TYPE_PROTEIN']
        filename = 'sars-cov2-genes.zip'
        api_response = api_instance.download_gene_package(gene_ids,
                            filename=filename, _preload_content=False)
        zip_file = open(filename, 'wb')
        try:
            copyfileobj(api_response, zip_file)
            with zipfile.ZipFile(filename) as zip_data:
                zip_data.extractall(dest)
        finally:
            zip_file.close()
            remove(filename)


def download_gene_data(gene_ids, dest):
    download_gene_data_via_api(gene_ids, dest)

def get_gene_data(gene_ids, dest):
    download_gene_data(gene_ids, dest)
    with open(f'{dest}/ncbi_dataset/data/data_report.yaml') as yaml_file:
        gene_data = yaml.load(yaml_file, Loader=yaml.SafeLoader)
    # rmtree(dest)
    return gene_data


def get_strand(orientation):
    """Convert NCBI Dataset data report orientation to BioPython strand.
    """
    return {'minus':-1, 'plus':+1, None:None}[orientation]


def make_simple_location(interval, ref, default_strand=None):
    """Convert NCBI Dataset data report range element to BioPython
    FeatureLocation.
    
    Note: Data report is 1-based, BioPython is 0-based.
    """
    strand = get_strand(interval.get('orientation', None))
    if strand is None:
        strand = default_strand
    return FeatureLocation(int(interval['begin'])-1,
                           int(interval['end']),
                           strand,
                           ref)


def make_location(location_data, default_strand=None):
    """Convert NCBI Dataset data report location to BioPython
    FeatureLocation or CompoundLocation.
    
    Note: Data report is 1-based, BioPython is 0-based.
    """

    ref = location_data['accessionVersion']
    if len(location_data['range']) == 1:
        return make_simple_location(location_data['range'][0],
                                    ref,
                                    default_strand)

    loc = []
    for interval in location_data['range']:
        loc.append(make_simple_location(interval, ref, default_strand))
    return CompoundLocation(loc)


def clip_location(position):
    """Clip a postion for a BioPython FeatureLocation to the reference.
    
    WARNING: Use of reference length not implemented yet.
    """
    return max(0, position)


def gene_to_upstream(gene, upstream_length = 2000):
    """Given a gene structure from an NCBI Dataset gene data report,
    return a list of rows with upstream range and the set of
    transcripts sharing that upstream range.
    
    Example: [[FeatureLocation(...), 'NM_123', 'NM_456']]
    """
    upstream_collection = dict()
    # Need strand from gene location, since it's not included in exon locations.
    gene_location = make_location(gene['genomicRanges'][0])

    for transcript in gene['transcripts']:
        transcript_accession = transcript['accessionVersion']
        transcript_location = make_location(transcript['exons'], gene_location.strand)

        if transcript_location.strand < 0:
            upstream = FeatureLocation(transcript_location.end,
                                       transcript_location.end + upstream_length,
                                       transcript_location.strand,
                                       gene_location.ref)
        else:
            upstream = FeatureLocation(clip_location(transcript_location.start - upstream_length),
                                       clip_location(transcript_location.start),
                                       transcript_location.strand,
                                       gene_location.ref)

        # BioSeq FeatureLocation not hashable.
        upstream_key = (min(upstream), max(upstream), upstream.strand, upstream.ref)
        if upstream_key in upstream_collection:
            upstream_collection[upstream_key].append(transcript_accession)
        else:
            upstream_collection[upstream_key] = [upstream, transcript_accession]
    return list(upstream_collection.values())


def write_bed(output, location, name, score=None):
    bed_strand = {None:'.', 0:'.', -1:'-', 1:'+'}[location.strand]
    bed_score = '.' if score is None else score
    output.write(f'{location.ref}\t{location.start}\t{location.end}\t{name}\t{bed_score}\t{bed_strand}\n')
