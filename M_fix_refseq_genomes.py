# Add strain name to entries in refseq multifastas, and replace spaces with underscores

import os
import fileinput
import csv

headers = open('bacteria_headers.txt','r')
header_cols = csv.reader(headers,delimiter='|')

unique_species = set()
for line in header_cols:
	name = line[4].lstrip('_').split('_')
	if 'sp.' in name:
		unique_species.add(tuple(name[0:3]))
	else:
		unique_species.add(tuple(name[0:2]))

print unique_species
print len(unique_species)
