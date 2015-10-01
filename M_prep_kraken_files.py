# Remove text in front of gi| in fasta header so kraken and clark can parse them.

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fasta_filename')
args = parser.parse_args()

mfa = open(args.fasta_filename,'r')
new_mfa = open(args.fasta_filename.rpartition('.')[0] + '_kraken.fasta','w')
for line in mfa:
	if line[0] == '>': # fasta name
		new_mfa.write('>gi|' + line.partition('gi|')[2]) # remove the strain label at the front of the header because kraken is stupid
	else:
		new_mfa.write(line)

new_mfa.close()
