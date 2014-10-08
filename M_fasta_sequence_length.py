#!/usr/bin/python
from Bio import SeqIO
import sys
import pprint
import csv

# Modified from http://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/

# Get fasta lengths
fasta_sizes = dict()
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
	fasta_sizes[seq_record.id] = len(seq_record)

# Format of sam header: @SQ, SN:<name>, LN:<length>
sam_header = csv.reader(open('testE_head.sam','r'),delimiter='\t')
for row in sam_header:
	if row[1][3:] not in fasta_sizes.keys():
		print row
	else:
		if int(row[2][3:]) != fasta_sizes[row[1][3:]]:

			print row
			print fasta_sizes[row[1][3:]]

