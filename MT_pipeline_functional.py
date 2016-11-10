# pull out fasta entries that match target_ids from kallisto output file and make new multifasta
import argparse
import csv

parser = argparse.ArgumentParser(description='Make mfa from target_ids from kallisto output')
parser.add_argument('filename', help='kallisto abundance.tsv file')
args = parser.parse_args()

# get gene IDs from kallisto output
with open(args.filename,'r') as input_file:
	input_csv = csv.reader(input_file, 'excel-tab')
	input_data = [r for r in input_csv]
	targets = set(g.partition('_')[0] for g in zip(*input_data)[0])

mfa = open('refined_mfa.fa','w')

with open('meta.nuc','r') as fasta:
	for line in fasta:
		if line.startswith('>'): #header line
			present = False
			if line.partition(' ')[0].partition('_')[0].partition('>')[2] in targets:
				present = True
				mfa.write(line)
		elif present:
			mfa.write(line)

mfa.close()
