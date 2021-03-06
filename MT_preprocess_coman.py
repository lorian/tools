# Preprocess transcript file to add strain names
import argparse
import csv
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'


mfa = open('processed_ReadCount_dataset2.txt','w')

with open('gene_ReadCount_dataset2.txt','r') as transcripts:
	for line in transcripts:
		if line.startswith('gi|'):
			# look up gi number: gi|479337003|ref|YP_007848221.1|
			uid = line.partition('gi|')[2].partition('|')[0]

			try:
				handle = Entrez.efetch("nucleotide", id=uid, retmode="xml")
				records = Entrez.read(handle)
			except:
				print "Error with {}".format(line)
				mfa.write(line)
				continue

			# standard location
			try:
				taxid = records[0]['GBSeq_feature-table'][0]['GBFeature_quals'][-1]['GBQualifier_value'].partition(':')[2]
			except:
				pass

			if not taxid:
				# sometimes it's not the last quals list, so we have to search for it
				try:
					taxid = [r['GBQualifier_value'] for r in records[0]['GBSeq_feature-table'][0]['GBFeature_quals'] if ('taxon' in r['GBQualifier_value'])][0].partition(':')[2]
				except:
					pass

			if taxid.isdigit():
				mfa.write('taxid|'+ taxid +'|'+ line)
			else:
				print "Error with {}".format(line)
				mfa.write(line)

		else:
			mfa.write(line)

mfa.close()
