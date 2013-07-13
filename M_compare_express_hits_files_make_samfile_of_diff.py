# Create sam file listing reads that are present in one express hits file but not another

import sys
import csv

filenameA = sys.argv[1]
filenameB = sys.argv[2]

fileA = open(filenameA,'r')
fileB = open(filenameB,'r')

namesA = set()
namesB = set()
for lineA in fileA:
	if not lineA.startswith('@'):
		aligned = lineA.split('\t')[2]
		if aligned != '*': #only keep aligned reads
			namesA.add(lineA[:lineA.find('\t')])

for lineB in fileB:
	if not lineB.startswith('@'):
		aligned = lineB.split('\t')[2]
		if aligned != '*': #only keep aligned reads
			namesB.add(lineB[:lineB.find('\t')])

fileA.close()
fileB.close()

diff = namesA ^ namesB
print "Items in file A: {0}".format(len(namesA))
print "Items in file B: {0}".format(len(namesB))
print "Items only in one file: {0}".format(len(diff))
print "Double-checking math: {0} + {1} = {2} = {3}".format(min(len(namesA),len(namesB)),len(diff),min(len(namesA),len(namesB))+len(diff),max(len(namesA),len(namesB)))

# Create file containing missing reads
if len(namesA) > len(namesB):
	bigfilename = filenameA
else:
	bigfilename = filenameB

bigfile = open(bigfilename,'r')
missing_reads = open('missingreads.sam','w')
missing_reads.seek(0) #move to beginning to overwrite current text

for line in bigfile:
	if line.startswith('@'):
		missing_reads.write(line)
	else:
		name = line[:line.find('\t')]
		if name in diff:
			missing_reads.write(line)
#			diff.remove(name) # don't include duplicate reads

missing_reads.truncate()
missing_reads.close()
bigfile.close()
