# Concatenate the entries of a multi-entry fasta file into one single entry, for all mfas in a directory

import sys
import os

dirname = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if dirname == "":
		dirname = arg #avoid extra space at beginning
	else:
		dirname = dirname + " " + arg

file_list = [f for f in os.listdir(dirname) if f.endswith('dna.genome.fa')]

for f in file_list:
	basename = f.partition('.dna.genome.fa')[0]
	mfa = open(os.path.join(dirname,f),'r')

	text = ""

	firstline = True
	for line in mfa:
		if firstline:
			if line.find('|') == -1:
				# add name to beginning of header
				text = '>' + basename.partition('.')[0].replace (" ", "_") + "|" + line[1:]
			else:
				text = line
			firstline = False
		elif line.startswith('>'):
			text += 'NNNNNNNNNN' #indicate possible gaps between chr, plasmids, shotgun pieces, etc
		elif line!= "\n": #skip empty lines within fasta
			text += line

	text = text + "\n" #add newline to end of file for eventual cat
	fa = open(os.path.join(dirname,basename + '.cat.fa'), 'w')
	fa.seek(0)
	fa.write(text)
	fa.truncate()
	fa.close()
