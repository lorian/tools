"""
Add taxonomy IDs to list of files excluded by CLARK
* Run M_prep_clark_files.py to create directory structure
* Copy printed output of file names
* Save as files_excluded.txt (or use the one created by CLARK when it runs)
* Copy directory names in Custom
* Replace underscores with spaces
* Paste into http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi to get tax_ids
* (Check "full taxid lineage")
* Save as file tax_report.txt
* Manually look up genuses for species that weren't auto-identified, save as file, and copy into tax_report.txt
"""

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
with open(os.path.expanduser('tax_report.txt'),'r') as tax:
	tax_csv = csv.reader(tax, 'excel-tab')
	next(tax_csv) # skip header line
	for line in tax_csv:
		if line[4] != ' ':
			if line[6].find(' ') != -1: # full taxid identification was successful
				tax_dict[lanthpy.single_name_cleanup(line[2])] = line[6].partition(" ")[2].partition(" ")[0] #[name] = tax_id
			else: # species only looked up
				tax_dict[lanthpy.single_name_cleanup(line[2])] = line[6]

pprint.pprint(tax_dict)

output = open(os.path.expanduser('i100_tax_ids.txt'),'w')

with open(os.path.expanduser('files_excluded.txt'),'r') as files:
	for line in files:
		tax_id = lookup_name(tax_dict,lanthpy.single_name_cleanup(line))
		if tax_id != 0:
			#print "{0}\t{1}\n".format(line.rstrip(),tax_id)
			output.write("{0}\t{1}\n".format(line.rstrip(),tax_id))
		else:
			print "Unable to lookup taxid for {}".format(lanthpy.single_name_cleanup(line))

