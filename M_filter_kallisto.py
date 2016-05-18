# go through read files, keep only reads that don't match pseudobam plus their next line

import csv
import argparse
import pprint
import collections

parser = argparse.ArgumentParser(description='Process kallisto pseudosam and remove all reads that pseudoaligned')
parser.add_argument('samname', help='kallisto pseudosam')
parser.add_argument('readsname', help='fastq containing original reads')
args = parser.parse_args()

with open(args.samname,'r') as sam_file:
	sam_csv = csv.reader(sam_file, delimiter='\t')
	sam_data = [r for r in sam_csv]
	matched_reads = set()
	for r in sam_data:
		if not r[0].startswith('@SQ'):
			matched_reads.add(r[0])

with open(args.readsname,'r') as read_file:
	read_csv = csv.reader(read_file, delimiter=' ')
	read_data = [r for r in read_csv]
	ready_to_delete = False
	for r in read_data:
		if len(r) > 1 and r[0][1:] in matched_reads: # remove @ from start of line
			print r[0]

		#if not r[0].startswith('@SQ'):
		#	matched_reads.add(r[0])


print len(matched_reads)
