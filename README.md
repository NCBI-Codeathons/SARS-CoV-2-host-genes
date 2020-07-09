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

- ACE2
- ABO
- TMPRSS2

Organism/gene table:

|organism           | tax_id | ACE2     | ABO      | TMPRSS2  |
|:----------------- |:------ |:-------- |:-------- |:-------- |
|Human              |9606    |59272     |28        |7113      |
|Mouse              |10090   |70008     |80908     |50528     |
|Zebrafish          |7955    |492331    |-         |494080    |
|Rhesus monkey      |9544    |712790    |722252    |715138    |
|Common vampire bat |9430    |112313373 |112320051 |112306012 |

## Sketch of output and how to produce it

Output | Description | How to present it
--- | --- | ---
location|genomic coordinates from most recent NCBI annotation|BED file
CDS|coordinates on the genome, and on the transcript|BED file
exons|coordinates on the genome|BED file
transcript variants|coordinates on the genome (start to stop), with transcript identifier|BED file
protein sequences|proteins from Protein, with NP identifier|FASTA
protein domain annotation|domains or protein subparts (e.g signal, mature peptides) as annotated on protein flat files; project to the genome to output with genomic coordinates too|FASTA for peptides, BED file for genomic coordinates
upstream regions|genomic coordinates of region 2 kb upstream of annotated transcription starts|BED file
UTRs|genomic coordinates for exonic regions of transcripts just before (5' UTR) and after (3' UTR) the CDS|BED file
introns|genomic coordinates for regions of transcript span that are not exonic|BED file
**Non-sequence metadata**|
summary|summary in Gene record|TSV file
gene symbols, aliases|primary gene symbol and aliases from Gene record|TSV file
primary description and other names|primary description and other names from Gene record|TSV file
publications|associated PubMed IDs from Gene record|TSV file
gene type|the type of gene from Gene record, e.g. protein-coding|TSV file
expression data|cell/tissue type expression data from Gene record|TSV file

## Technical Implementation
- Prototyping of data retrieval
    - Command-line `datasets` tool
    - `jq`
    - Simple shell scripting
- Python for data retrieval and transformation
    - Datasets API for data retrieval
    - Use a library for generating BED files
    - A top-level Python script for coordinating everything
    - User can input a list of GeneIDs to run the script
- Visualize via web page or Jupyter Notebook
    - Aiming to reuse/embed GDV to show multiple genomes, aligned, on one page
    
### BED file specs
See http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for list of standard columns and specifications. Here we use:

- 6-column BED format given the need to capture strand for genes and subfeatures
- RefSeq accession.version for column 1; also more explicit when multiple organisms are included in output
- 0-based start coordinates, 1-based end coordinates
- File sorting by column 1 primary, column 2 secondary (chr and start positions, `sort -k1,1 -k2,2n` on command line)
- Score (column 5) is irrelevant to this output but requires an integer, so value 1 is used here
- The GeneID is indicated in column 4 (the only descriptive/label column), and may be joined via an underscore to other pertinent identifiable information for the feature. The name/label does not include spaces.
- The relevant strand in column 6 is indicated by either '+' or '-'. All gene-related features in this output are strand-specific.

### Non-sequence metadata
* TSV file with one gene per row
* Extracted from public resouces: esummary, efetch, elink
* Progress:
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

### Datasets tool
See https://www.ncbi.nlm.nih.gov/datasets/docs/command-line-start/ for documentation on how to download and use this tool.

To obtain the specific GeneIDs used in this Codeathon (as in organism/gene table above), use:

```./datasets gene-descriptors gene-id 59272 28 7113 70008 80908 50528 492331 494080 712790 722252 715138 112313373 112320051 112306012|jq```

To download these genes to a file use:

```./datasets download gene 59272 28 7113 70008 80908 50528 492331 494080 712790 722252 715138 112313373 112320051 112306012```

### Sample output
... to be continued.

