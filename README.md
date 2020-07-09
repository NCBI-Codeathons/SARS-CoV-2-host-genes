# Characterization of SARS-CoV-2 host genes
Datasets Codeathon Team II.

### Team members
Andrew Kim  
Bart Trawick  
Catherine Farrell  
Pooja Strope  
Wratko Hlavina  
Xuan Zhang  

## Introduction
The study of SARS-CoV-2 has become a significant interest for human health research, and has drawn attention to the human genes associated with SARS-CoV-2 entry. Users are therefore interested in retrieving as much information as possible about these host genes and their products, e.g. ACE2 encoding the SARS-CoV-2 spike protein receptor or TMPRSS2 encoding a protease that facilitates viral entry. In order to facilitate study of SARS-CoV-2 host genes, a comprehensive report of key features of structural annotation will be produced from this Codeathon, and will also include key metadata for these genes, such as nomenclature, summaries or associated publications. Since study of these genes in model organisms is also of interest to researchers, information on the orthologous genes in a select set of organisms will also be produced. The overriding goal will be to provide data from various NCBI resources (e.g. Gene, Nucleotide, Protein records) in a succinct tabular format, such that the user can quickly retrieve multiple types of high-interest information about these genes in one place.

## Scope
Genes to work on:

- ACE2 GeneID:59272
- ABO GeneID:28
- TMPRSS2 GeneID:7113

## Progress

#### Non-sequence metadata
* Extracted from esummary, efetch, elink
* One gene
* Sample output:
```angular2
Summary: This gene encodes proteins related to the first discovered blood group system, ABO. Variation in the ABO gene (chromosome 9q34.2) is the basis of the ABO blood group, thus the presence of an allele determines the blood group in an individual. The 'O' blood group is caused by a deletion of guanine-258 near the N-terminus of the protein which results in a frameshift and translation of an almost entirely different protein. Individuals with the A, B, and AB alleles express glycosyltransferase activities that convert the H antigen into the A or B antigen. Other minor alleles have been found for this gene. [provided by RefSeq, Sep 2019]
Symbol: ABO
Aliases: A3GALNT, A3GALT1, GTB, NAGAT
Description: ABO, alpha 1-3-N-acetylgalactosaminyltransferase and alpha 1-3-galactosyltransferase
Type:  protein-coding
Publications: 32093636, 31595994, 31487040, 31329303, 31260107, 31240718, 30859643, 30791881, 30549285, 30347622, 29873711, 29659952, 29457687, 29106707, 28984382, 28939368, 28833251, 28653406, 28256109, 27979997, 27542834, 27538125, 26924317, 26632894, 26512559, 26329815, 26268879, 26247473, 26244499, 26148378, 25854361, 25820620, 25656610, 25636112, 25297604, 25217989, 25156869, 25138306, 25064734, 24743543, 24510570, 23816557, 23387832, 23319424, 23300549, 22963146, 22642827, 22258027, 21729554, 21560847, 21306478, 20703243, 20677133, 20666915, 20576794, 20456702, 20371059, 20197725, 20154292, 20003128, 19648918, 19490215, 19470260, 19276450, 19078892, 19054377, 18712158, 18680548, 18651204, 18426679, 18273824, 18247104, 18156754, 18078207, 18067076, 18063521, 18003641, 17642512, 17531777, 17393014, 17311872, 17259183, 17130965, 17002642, 16871363, 16686846, 16686845, 16631357, 16403294, 16239542, 16215642, 16181218
```