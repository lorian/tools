# Add taxonomy IDs to list of files excluded by CLARK

import csv
import pprint
import os
import lanthpy

def lookup_name(sp_dict,line):
	for sp in sp_dict.keys():
		if sp in line:
			#print "{0}: {1} -> {2}".format(line,sp,sp_dict[sp])
			return sp_dict[sp]
	return 0

tax_dict = dict()
with open(os.path.expanduser('~/metagenome/i100_species_tax.txt'),'r') as tax:
	tax_csv = csv.reader(tax, 'excel-tab')
	for line in tax_csv:
		tax_dict[lanthpy.single_name_cleanup(line[2])] = line[6]

pprint.pprint(tax_dict)

output = open(os.path.expanduser('~/metagenome/i100_tax_ids.txt'),'w')
with open(os.path.expanduser('~/metagenome/files_excluded.txt'),'r') as files:
	for line in files:
		tax_id = lookup_name(tax_dict,lanthpy.single_name_cleanup(line))
		if tax_id != 0:
			output.write("{0}\t{1}\n".format(line.rstrip(),tax_id))
