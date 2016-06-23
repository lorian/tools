# go through read files, keep only reads that match a given pattern
# Note: read files must be checked individually

import csv
import argparse
import pprint
import collections
import os

parser = argparse.ArgumentParser(description='Create fastq containing only specific reads')
parser.add_argument('keepreads', help='Phrase fastq header must contain in order to be kept')
parser.add_argument('readsname', help='fastq containing original reads')
args = parser.parse_args()

new_reads = open(args.readsname.partition(".")[0] + "_{}.fastq".format(keepreads),"w")

with open(args.readsname,'r') as read_file:
	read_csv = csv.reader(read_file, delimiter=' ')
	read_data = [r for r in read_csv]
	ready_to_keep = False
	for r in read_data:
		if len(r) > 1: #read header
			if keepreads in r[0][1:]: # remove @ from start of line
				ready_to_keep = True
			else:
				ready_to_keep = False

		if ready_to_keep:
			new_reads.write(" ".join(r,).strip()+"\n")

new_reads.close()
print len(matched_reads)


