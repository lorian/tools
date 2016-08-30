"""
Add taxonomy IDs to list of files excluded by CLARK
* Run from clark db directory
* Get files_excluded.txt created by CLARK when it runs
* Get taxids from headers in files in Custom
* Manually look up genuses for species that weren't auto-identified, save as file, and copy into tax_report.txt
"""

import csv
import pprint
import os
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'

def get_taxid(f):
		# get taxid
	mfa = open(os.path.join('Custom/',f),'r')
	
	for line in mfa:
		taxid = line.partition('taxid|')[2].partition('|')[0]
		return taxid
		
original_files = [f for f in os.listdir('Custom/') if f.endswith(".cat.fa")]

taxids = dict()
for f in original_files:
	taxid = get_taxid(f)
	taxids[f] = taxid
	'''
	# lookup full tax lineage
	handle = Entrez.efetch(db="taxonomy", id=taxid, mode="text", rettype="xml")
	taxon = Entrez.read(handle)[0] # only grab the first

	print '\n'
	print taxon["LineageEx"]
	
	lineage = []
	try:
		for t in taxon["LineageEx"]:
			lineage.append(t['TaxId'])
		if taxid not in lineage:
			lineage.append(taxid)
	except:
		pass # no parent taxonomy

	print lineage
	'''
	
output = open(os.path.expanduser('targets.txt'),'w')

with open(os.path.expanduser('files_excluded.txt'),'r') as files:
	next(files) # skip first line
	for line in files:
		taxid = taxids[line.rpartition('/')[2].rstrip()]
		if taxid:
			print "{0}\t{1}\n".format(line.rstrip(),taxid)
			output.write("{0}\t{1}\n".format(line.rstrip(),taxid))
		else:
			print "Unable to lookup taxid for {}".format(line)
