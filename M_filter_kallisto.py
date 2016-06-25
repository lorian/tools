# go through read files, keep only reads that don't match pseudobam plus their next line
# Note: read files must be checked individually

import csv
import argparse
import pprint
import collections
import os

parser = argparse.ArgumentParser(description='Process kallisto pseudosam and remove all reads that pseudoaligned')
parser.add_argument('samname', help='kallisto pseudosam')
parser.add_argument('readsname', help='fastq containing original reads')
args = parser.parse_args()

'''
mapped_sam = args.samname.partition('.')[0] + "_mapped.sam"
if not os.path.isfile(mapped_sam):
	# Only keep mapped reads
	os.system('samtools view -S -F 4 {} > {}'.format(args.samname,mapped_sam))
'''
# use raw sam:
mapped_sam = args.samname

with open(mapped_sam,'r') as sam_file:
	sam_csv = csv.reader(sam_file, delimiter='\t')
	sam_data = [r for r in sam_csv]
	matched_reads = set()
	for r in sam_data:
		if not r[0].startswith('@SQ'):
			matched_reads.add(r[0])

new_reads = open(args.readsname.partition(".")[0] + "_filtered_raw.fastq","w")

with open(args.readsname,'r') as read_file:
	read_csv = csv.reader(read_file, delimiter=' ')
	read_data = [r for r in read_csv]
	ready_to_delete = False
	for r in read_data:
		if len(r) > 1: #read header
			if r[0][1:] in matched_reads: # remove @ from start of line
				ready_to_delete = True
			else:
				ready_to_delete = False

		if not ready_to_delete:
			new_reads.write(" ".join(r,).strip()+"\n")

new_reads.close()
print len(matched_reads)


