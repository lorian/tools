#!/usr/bin/python
from Bio import SeqIO
import sys

# Modified from http://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/

cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
	output_line = "{0}\t{1}".format(seq_record.id, len(seq_record))
	print(output_line)
