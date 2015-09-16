# Counts length of fastas in a multi-fasta file
# from https://bioexpressblog.wordpress.com/2014/04/15/calculate-length-of-all-sequences-in-an-multi-fasta-file/

from Bio import SeqIO
import sys

cmdargs = str(sys.argv)
fasta_name = str(sys.argv[1])
output_file = open('{}_lengths.txt'.format(fasta_name.rpartition('.')[0]),'w')
for seq_record in SeqIO.parse(fasta_name, "fasta"):
	output_line = '%s\t%i' % (seq_record.id, len(seq_record))
	output_file.write(output_line+'\n')
