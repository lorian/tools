# Print all multimapped reads that align with indels

import sys
import csv
import pysam

filename = sys.argv[1]
sam = open(filename,'r')

last_read = ""
last_line = ""
contents = ""
for line in sam:
	if not line.startswith('@'):
		fields = line.split('\t')
		qname = fields[0]
		if qname == last_read: #multi-mapped read
			cigar = fields[5]
			if 'S' in cigar or 'I' in cigar or 'H' in cigar or 'D' in cigar:
#				print "{0} aligns as {1}".format(qname,cigar)
				contents = contents + line
		last_read = qname
		last_line = line
	if last_read.startswith('B'):
		break

sam.close()

multireads = open('multi_indel_reads.txt','w')
multireads.seek(0) #move to beginning to overwrite current text
multireads.write(contents)
multireads.truncate()
multireads.close()

