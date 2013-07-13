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

file_list = os.listdir(dirname)

for f in file_list:
	if f.find('.mfa') != -1:
		basename = f[:-4]
		mfa = open(os.path.join(dirname,f),'r')

		text = ""

		firstline = 1
		for line in mfa:
			if firstline == 1:
#				text = '>' + basename.replace (" ", "_") + " " + line[1:] #add name to beginning of ID line
				text = line
				firstline = 0
			elif line.startswith('>'):
				text = text + '-'
			elif line!= "\n": #skip empty lines within fasta
				text = text + line

		text = text + "\n" #add newline to end of file for eventual cat
		fa = open(os.path.join(dirname,basename + '.fa'), 'w')
		fa.seek(0)
		fa.write(text)
		fa.truncate()
		fa.close()
