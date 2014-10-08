#!/usr/bin/python
from Bio import SeqIO
import sys
import pprint

# Modified from http://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/

fasta_sizes = dict()
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
	fasta_sizes[seq_record.id] = len(seq_record)

print fasta_sizes['EUKY_12|gi|134108779|ref|NC_009178.1|_Cryptococcus_neoformans_var._neoformans_B-3501A_chromosome_2,_whole_genome_shotgun_sequence913558']
