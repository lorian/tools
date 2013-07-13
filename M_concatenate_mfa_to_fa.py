# Concatenate the entries of a multi-entry fasta file into one single entry

import sys

filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if filename == "":
		filename = arg #avoid extra space at beginning
	else:
		filename = filename + " " + arg

basename = filename[:filename.rfind('.')]

mfa = open(filename,'r')
text = ""

triggered = 0
firstline = 1
for line in mfa:
	if firstline == 1:
#		text = '>' + basename.replace (" ", "_") + " " + line[1:] #add name to beginning of ID line
		text = line
		firstline = 0
	elif triggered == 1:
		if line.startswith('>'):
			text = text + '-'
		triggered = 0
	elif line == "\n":
		triggered = 1
	else:
		text = text + line

fa = open(filename[:filename.rfind('.')] + '.fa', 'w')
fa.seek(0)
fa.write(text)
fa.truncate()
fa.close()
