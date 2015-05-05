# Look through a directory of mfa's for those that lack an actual genome sequence
import sys
import os

def has_target(entries, target): #is there target phrase in the list?
	for e in entries:
		if e.find('plasmid') == -1 and e.find(target) != -1:
			return True
	return False


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
	if f.endswith('.mfa'):
		basename = f[:-4]
		mfa = open(os.path.join(dirname,f),'r')
		entries = []

		for line in mfa:
			if line.startswith('>'):
				entries.append(line)

		if not has_target(entries, 'complete genome') and not (len(entries) == 1 and has_target(entries,'chromosome')):
			print "{0} lacks genome".format(basename)
			for e in entries:
				print "\t{0}".format(e)

