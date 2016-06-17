'''
Mash results are tab delimited lists of Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes:
Bacillus_aryabhattai_gca_001043825.GCA_001043825.1.30.dna.genome.fa     /home/lorian/scratch/illumina_100species_trimmed.1.fq   1       1       0/1000
Bacillus_atrophaeus_1942.GCA_000165925.1.30.dna.genome.fa       /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.295981        0.000944763     1/1000
Bacillus_bombysepticus_str_wang.GCA_000831065.1.30.dna.genome.fa        /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.197292        2.3647e-28      8/1000
'''

import csv
import argparse
import pprint
import collections
import os

parser = argparse.ArgumentParser()
parser.add_argument('filename', help='mash output file')
args = parser.parse_args()

with open(args.filename,'r') as mash_file:
	mash_csv = csv.reader(mash_file, delimiter='\t')
	mash_data = [r for r in mash_csv]
	for r in mash_data:
		if r[0] != r[1] and r[2] != '1': # ignore self-self comparisons, and completely unrelated genomes
			print r
