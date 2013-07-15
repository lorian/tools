# Process large metagenomic reference dataset:
#	Make names fully readable by replacing spaces with _
#	Discard virus genomes
#	Combine all entries for a single species using 10 N's
#	Split into files of less than 3.6 billion characters for bowtie2 index

import sys

filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if filename == "":
		filename = arg #avoid extra space at beginning
	else:
		filename = filename + " " + arg

mfa = open(filename,'r')
text = ""

for line in mfa:
	if line.startswith('>'):
		text = text + line.replace (" ", "_")
	else:
		text = text + line

fa = open(filename[:filename.rfind('.')] + '_processed.fa', 'w')
fa.seek(0)
fa.write(text)
fa.truncate()
fa.close()
