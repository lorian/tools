# Run against mash index in parallel

import subprocess
import string
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fastq', help='Fastq file to calculate distance from mash indexes')
args = parser.parse_args()

for letter in string.ascii_uppercase+"_":
	print letter
	subprocess.call('time mash dist ~/scratch/ensembl_bact/ignore/ensembl_bact_{0}.msh {2} > mash_real_{1}_{0}.txt &'.format(letter, args.fastq.partition('.')[0], args.fastq), shell=True)

subprocess.call('cat mash_real_{0}_* > mash_real_{0}.txt'.format(args.fastq.partition('.')[0]), shell=True)
