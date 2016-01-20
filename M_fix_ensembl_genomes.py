# Add strain name to entries in ensembl multifasta, and replace spaces with underscores

import os
import fileinput

dirname = '.'
fasta_list = [f for f in os.listdir(dirname) if f.endswith('.fa')]

for f in fasta_list:
	name = f.partition('.')[0] # get strain name from filename
	print name
	for line in fileinput.input(os.path.join(dirname,f), inplace=True):
		if line.startswith('>') and not line.startswith('>'+ name +'|'):
			print '>'+ name + "|" + line[1:].replace(" ","_"), # writes to file: note comma!
		else:
			print line,

