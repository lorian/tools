# Module containing frequently-used functions

import csv
import string
import os
import sys

def get_script_path():
	return os.path.dirname(os.path.realpath(sys.argv[0]))

def genome_name_cleanup(raw_names):
	"""
	Takes a list of messy genome names such as found in the Martin_etal
	multifasta and returns cleaned-up, normalized versions.
	"""

	#if BACT_ abbreviations are found, see if they match large genome database
	try:
		min( i for i, sp in enumerate(raw_names) if 'BACT_' in sp )
	except:
		clean_names = raw_names
	else:
		with open(os.path.join(get_script_path(),'M_biggenomenames.txt'), 'r') \
					as biggenomefile:
			big_genome = dict(csv.reader(biggenomefile))
		clean_names = [big_genome[s.partition('|')[0]]
						if s.partition('|')[0] in big_genome.keys()
						else s for s in raw_names]

	# Force clean_names to be lowercase, underscored, and without punctuation
	bad_punct = "!\"#$%&'()*+,-./:;<=>?@[\]^`{|}~" # string.punctuation w/no _
	clean_names = [s.lower().translate(string.maketrans("",""),bad_punct)
					.replace(" ","_") for s in clean_names]

	return clean_names
