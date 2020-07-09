# This script will print "cds" from datasets download

import yaml
import pprint

pp = pprint.PrettyPrinter(indent=4)

a_yaml_file = open("data_report.yaml")
bed_file = open('cds_output.bed', 'a')
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

cds = parsed_yaml_file["genes"][0]["transcripts"][0]["cds"]
chromosome = parsed_yaml_file["genes"][0]['chromosomes'][0]
accessionVersion = parsed_yaml_file["genes"][0]['genomicRanges'][0]['accessionVersion']
chromStart = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["begin"]
chromEnd = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["end"]
orientation = parsed_yaml_file["genes"][0]["genomicRanges"][0]["range"][0]["orientation"]

pp.pprint(cds)
# prints out the following...
# {'accessionVersion': 'NM_021804.3', 'range': [{'begin': '307', 'end': '2724'}]}

bed_file.write(f"{accessionVersion}    {chromStart}    {chromEnd}\n")
a_yaml_file.close()
bed_file.close()