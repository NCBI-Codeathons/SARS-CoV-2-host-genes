# This script will print "cds" from datasets download

import yaml
import pprint

pp = pprint.PrettyPrinter(indent=4)

a_yaml_file = open("data_report.yaml")
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

pp.pprint(parsed_yaml_file["genes"][0]["transcripts"][0]["cds"])

# prints out the following...
# {'accessionVersion': 'NM_021804.3', 'range': [{'begin': '307', 'end': '2724'}]}