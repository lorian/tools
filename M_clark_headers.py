# run from inside directory with genomes to edit header of; will add gi and gca to front, if they exist.

import csv
import argparse
import pprint
import collections
import os
import re

def change_header(f):
	basename = f.partition('.cat.fa')[0]
	try:
		mfa = open(f,'r')
	except IOError:
		print "{} does not exist".format(f)
		return False

	text = ""

	firstline = True
	for line in mfa:
		if firstline:
			firstline = False
			# replace header with name
			new_name = basename.partition('.1.30')[0].partition('.GCA')[0].replace(" ", "_")
			# try to filter out only the useful parts of the original header: GCA, NC, or gi uids
			name_parts = filter(None, re.split("[|:]+", line.replace('gi|','gi_').replace('_gi', '|gi')))

			# gi|XXXXXX|gb|ACCESSION.VERSION| Name for CLARK
			NC = [w for w in name_parts if (w.startswith('NC_') or w.startwith('AC_'))]
			GI = [w for w in name_parts if w.startswith('gi_')]
			if NC and GI:
				text = ">gi|{0}|gb|{1}|{2}\n".format(NC[0],GI[0],line[1:])
			elif NC:
				text = ">{0}|{1}\n".format(NC[0],line[1:])
			elif GI:
				text = ">gi|{0}|{1}\n".format(GI[0],line[1:])
			else:
				text = line
				print "Couldn't get NC for {}".format(f)
				print line
	
		elif line.startswith('>'):
			text += 'NNNNNNNNNN' #indicate possible gaps between chr, plasmids, shotgun pieces, etc
		elif line!= "\n": #skip empty lines within fasta
			text += line

	text = text + "\n" #add newline to end of file for eventual cat
	fa = open(basename + '.clark.fa', 'w')
	fa.seek(0)
	fa.write(text)
	fa.truncate()
	fa.close()
	return basename + '.clark.fa'

# Run in directory containing files to edit headers
original_files = [f for f in os.listdir('.') if f.endswith(".cat.fa")]

for f in original_files:
	change_header(f)
	
	
