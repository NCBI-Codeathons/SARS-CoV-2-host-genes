# This script will print "exons" from datasets download

import yaml
import pprint

pp = pprint.PrettyPrinter(indent=4)

a_yaml_file = open("data_report.yaml")
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

pp.pprint(parsed_yaml_file["genes"][0]["transcripts"][0]["exons"])

# prints out the following...
# {   'accessionVersion': 'NC_000023.11',
#     'range': [   {'begin': '15601956', 'end': '15602158', 'order': 1},
#                  {'begin': '15600726', 'end': '15601014', 'order': 2},
#                  {'begin': '15594845', 'end': '15595003', 'order': 3},
#                  {'begin': '15592229', 'end': '15592322', 'order': 4},
#                  {'begin': '15591713', 'end': '15591856', 'order': 5},
#                  {'begin': '15589344', 'end': '15589456', 'order': 6},
#                  {'begin': '15587753', 'end': '15587858', 'order': 7},
#                  {'begin': '15585475', 'end': '15585572', 'order': 8},
#                  {'begin': '15581221', 'end': '15581390', 'order': 9},
#                  {'begin': '15578089', 'end': '15578315', 'order': 10},
#                  {'begin': '15575666', 'end': '15575810', 'order': 11},
#                  {'begin': '15573367', 'end': '15573465', 'order': 12},
#                  {'begin': '15572201', 'end': '15572323', 'order': 13},
#                  {'begin': '15571624', 'end': '15571796', 'order': 14},
#                  {'begin': '15570295', 'end': '15570353', 'order': 15},
#                  {'begin': '15567726', 'end': '15567826', 'order': 16},
#                  {'begin': '15566253', 'end': '15566369', 'order': 17},
#                  {'begin': '15564024', 'end': '15564218', 'order': 18},
#                  {'begin': '15561033', 'end': '15562013', 'order': 19}]}