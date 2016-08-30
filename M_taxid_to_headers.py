# Collapse mfa contigs, look up and add taxid to header info in the run directory
import collections
import os
import re
import urllib2
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'

def get_taxid(uid):
	# look up taxid
	if uid.startswith('GCA_'):
		page_taxa = urllib2.urlopen('http://www.ebi.ac.uk/ena/data/view/{}&display=xml'.format(uid))
		for line in page_taxa:
			if line.strip().startswith('<TAXON_ID>'):
				return line.partition('>')[2].partition('<')[0]

	if uid.startswith('NC_') or uid.startswith('gi_'):
		if uid.startswith('gi_'):
			uid = uid.partition('gi_')[2] #NCBI wants them as just straight numbers
		handle = Entrez.efetch("nucleotide", id=uid, retmode="xml")
		records = Entrez.read(handle)
		# standard location
		taxid = records[0]['GBSeq_feature-table'][0]['GBFeature_quals'][-1]['GBQualifier_value'].partition(':')[2]
		if not taxid:
			# sometimes it's not the last quals list, so we have to search for it
			taxid = [r['GBQualifier_value'] for r in records[0]['GBSeq_feature-table'][0]['GBFeature_quals'] if ('taxon' in r['GBQualifier_value'])][0].partition(':')[2]
		return taxid

	return ''

def collapse_contigs(f):
	basename = f.partition('.dna.genome.fa')[0].partition('.mfa')[0]
	try:
		mfa = open(f,'r')
	except IOError:
		print "{} does not exist".format(f)
		return False

	text = ""

	firstline = True
	for line in mfa:
		if firstline:
			# replace header with name
			new_name = basename.partition('.1.30')[0].partition('.GCA')[0].replace(" ", "_")
			# try to filter out only the useful parts of the original header: GCA, NC, or gi uids
			name_parts = filter(None, re.split("[|:]+", line.replace('gi|','gi_').replace('_gi', '|gi')))
			keep_ids = [w for w in name_parts if (w.startswith('GCA_') or w.startswith('NC_') or w.startswith('gi_'))]
			for id in keep_ids:
				taxid = get_taxid(id)
				if taxid:
					break
			if not taxid:
				print "Unable to lookup taxid for {}".format(f)
				print line
				print keep_ids
			
			text = ">{0}|kraken:taxid|{1}|{2}\n".format(new_name,taxid,"|".join(keep_ids))
			firstline = False
		elif line.startswith('>'):
			text += 'NNNNNNNNNN' #indicate possible gaps between chr, plasmids, shotgun pieces, etc
		elif line!= "\n": #skip empty lines within fasta
			text += line

	text = text + "\n" #add newline to end of file for eventual cat
	fa = open(basename + '.cat.fa', 'w')
	fa.seek(0)
	fa.write(text)
	fa.truncate()
	fa.close()
	return basename + '.cat.fa'

original_files = [f for f in os.listdir('.') if f.endswith("fa")]

for f in original_files
	collapse_contigs(f)
