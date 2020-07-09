import yaml
import pprint


def main():
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


if __name__ == "__main__":
    main()
