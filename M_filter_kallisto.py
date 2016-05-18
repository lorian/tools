# import pseudosam
# pull out read names
# go through read files, keep only reads that don't match plus their next line
#"grep -v \"{}\" samname.fasta".format("|".join(readnames))


import csv
import argparse
import pprint
import container

parser = argparse.ArgumentParser(description='Process kallisto pseudosam and remove all reads that pseudoaligned')
parser.add_argument('samname', help='kallisto pseudosam')
args = parser.parse_args()

with open(args.samname,'r') as sam_file:
	sam_csv = csv.reader(sam_file, delimiter='\t')
	sam_data = [r for r in sam_csv]
	matched_reads = []
	for r in sam_data:
		print r[0]
		#matched_reads.append()
