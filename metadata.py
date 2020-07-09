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

class GeneRetrieve():
    def __init__(self, gene_id):
        Entrez.email = "who@nih.gov"
        Entrez.api_key = "ab0568529a7dd0e599fd12b3498f1c8e9e08"
        self.handle_summary = Entrez.esummary(db="gene", id=gene_id, rettype="xml")
        self.handle_full = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
        self.handle_gene2pubmed = Entrez.elink(dbfrom="gene", db="pubmed", id=gene_id, rettype="xml")
        # webUrl = urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=28')
        # webUrl.getcode() == 200
        # data = webUrl.read()

    def record(self):
        record_summary = Entrez.read(self.handle_summary)
        record_full = Entrez.read(self.handle_full)
        record_gene2pubmed = Entrez.read(self.handle_gene2pubmed)
        return record_summary['DocumentSummarySet']['DocumentSummary'][0], record_full[0], record_gene2pubmed[0]


class GeneData():
    def __init__(self, gene_id):
        gene_record = GeneRetrieve(gene_id)
        self.doc_summary, self.doc_full, self.doc_pubmed = gene_record.record()
        # self.doc = record['DocumentSummarySet']['DocumentSummary'][0]

    def pretty_print_field(self, dom_element):
        if dom_element.in_esummary == "Gene2Pubmed":
            return dom_element.pretty_print(self.doc_pubmed)
        if dom_element.in_esummary:
            return dom_element.pretty_print(self.doc_summary)
        return dom_element.pretty_print(self.doc_full)


class DomElement():
    def __init__(self, tag_name, attribute_name=None,  in_esummary=True):
        self.tag = tag_name
        self.attribute_name = attribute_name
        self.in_esummary = in_esummary

    def pretty_print(self, doc):
        if self.in_esummary == "Gene2Pubmed":
            pubmed_ids = []
            for link in doc['LinkSetDb'][0]['Link']:
                pubmed_ids.append(link['Id'])
            return ', '.join(pubmed_ids)
        if self.attribute_name:
            return getattr(doc[self.tag], 'attributes')[self.attribute_name]
        else:
            return doc[self.tag]


class GeneMetaReport():

    def __init__(self, gene_list):
        self.genes = gene_list
        self.fout_name = "metadata.tsv"
        self.field_list = [('Summary', DomElement('Summary', None, True)),
                  ('Symbol', DomElement('Name', None, True)),
                  ('Aliases', DomElement('OtherAliases', None, True)),
                  ('Description', DomElement('Description', None, True)),
                  ('Type', DomElement('Entrezgene_type', 'value', False)),
                  ('Publications', DomElement('', None, "Gene2Pubmed")),
                 ]

    def write_report(self):
        with open(self.fout_name, "w") as fout:
            fields = []
            for field in self.field_list:
                fields.append(field[0])
            fout.write('\t'.join(fields) + "\r\n")
            for gene_id in self.genes:
                data = GeneData(gene_id)
                field_val = []
                for _, val in self.field_list:
                    field_val.append(data.pretty_print_field(val))
                fout.write('\t'.join(field_val) + "\r\n")
        fout.close()
