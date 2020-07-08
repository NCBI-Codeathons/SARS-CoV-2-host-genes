# This script will print "name" from datasets download

import yaml
import pprint

pp = pprint.PrettyPrinter(indent=4)

a_yaml_file = open("data_report.yaml")
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

pp.pprint(parsed_yaml_file["genes"][0]["transcripts"][0]["name"])

# prints out the following...
# 'transcript variant 2'