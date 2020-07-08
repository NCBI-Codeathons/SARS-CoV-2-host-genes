"""
This Python retrieves Non-sequence metadata
Input: gene id
Output: print stdout
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=28
"""
import xml.etree.ElementTree as ET
from xml.dom import minidom
from Bio import Entrez
import urllib.request

HTTP_PREFIX = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&retmode=xml&id="


class GeneFetch():
    def __init__(self, gene_id):
        Entrez.email = "codeathon@example.com"
        self.handle_summary = Entrez.esummary(db="gene", id=gene_id, rettype="xml")
        self.handle_full = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        # webUrl = urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=28')
        # webUrl.getcode() == 200
        # data = webUrl.read()

    def record(self):
        record_summary = Entrez.read(self.handle_summary)
        record_full = Entrez.read(self.handle_full)
        return record_summary['DocumentSummarySet']['DocumentSummary'][0], record_full[0]
        # dom = minidom.parseString(data)


class GeneData():
    def __init__(self, gene_id):
        gene_record = GeneFetch(gene_id)
        self.doc_summary, self.doc_full = gene_record.record()
        # self.doc = record['DocumentSummarySet']['DocumentSummary'][0]

    def print_field(self, name, mapped_key):
        in_summary, field_name = mapped_key
        if in_summary:
            print(name + '\t' + self.doc_summary[field_name])
        # e = dom.getElementsByTagName('Entrezgene_summary')
        # print(e[0].childNodes[0].data)


if __name__ == "__main__":
    data = GeneData(28)
    value_list = {'Summary': (True, 'Summary'),
                  'Symbol': (True, 'Name'),
                  'Aliases': (True, 'OtherAliases'),
                  'Description': (True, 'Description'),
                  }
    for key, val in value_list.items():
        data.print_field(key, val)
