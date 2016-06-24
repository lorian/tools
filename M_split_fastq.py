# Split i100 fastq into individual fastqs based on source of reads

import csv
import argparse
import pprint
import collections
import os

parser = argparse.ArgumentParser(description='Split fastq by read source')
parser.add_argument('readsname', help='fastq containing original reads')
args = parser.parse_args()

with open(args.readsname,'r') as read_file:
	read_csv = csv.reader(read_file, delimiter=' ')
	read_data = [r for r in read_csv]
	print "Processing data..."
	split_data = collections.defaultdict(list)
	current_source = ""
	for r in read_data:
		if r[0].startswith('@'): #read header
			current_source = r[0][1:].partition('.fna')[0]
			split_data[current_source].append(" ".join(r,).strip()+"\n")
		else:
			split_data[current_source].append(" ".join(r,).strip()+"\n")

pprint.pprint(split_data.keys())
for k in split_data.keys():
	new_reads = open(args.readsname.partition(".")[0] + "_{}.fastq".format(k),"w")
	new_reads.write("".join(split_data[k],))
	new_reads.close()


