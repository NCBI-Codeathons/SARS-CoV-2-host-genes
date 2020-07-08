"""
This Python retrieves Non-sequence metadata
Input: gene id
Output: print stdout
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=28
"""
import xml.etree.ElementTree as ET
from xml.dom import minidom
from Bio import Entrez

HTTP_PREFIX = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id="


class GeneSummary():
    def __init__(self, gene_id):
        Entrez.email = "codeathon@example.com"
        self.handle = Entrez.esummary(db="gene", id=2)
        # Entrez.esummary(db="gene", id="28")

    def record(self):
        record = Entrez.read(self.handle, validate=False)
        return record


class GeneData():
    def __init__(self, gene_id):
        gene_sum = GeneSummary(gene_id)
        record = gene_sum.record()
        self.doc = record['DocumentSummarySet']['DocumentSummary'][0]

    def summary(self):
        summary = self.doc['Summary']
        print(summary)

if __name__ == "__main__":
    data = GeneData(28)
    data.summary()
