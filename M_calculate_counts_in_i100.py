# Calculate counts for each species in true i100 dataset

import argparse
import csv
import os
import collections
import cPickle
import pprint

parser = argparse.ArgumentParser()
parser.add_argument('genome_filename') # file linking reads to genomes
parser.add_argument('fastq_filename')
args = parser.parse_args()

# Set up dictionary converting between IDs and genome names
genome_csv = csv.reader(open(args.genome_filename,'r'), 'excel')
genome_data = [r for r in genome_csv if (r[1].startswith('NC') or r[1].startswith('AC'))]
genome_ids = list(zip(*genome_data)[1])
genome_names = list(zip(*genome_data)[3])

genome_data_addl = [r for r in genome_data if r[2].startswith('NC')] # second column contains second chromosomes, occasionally
genome_ids.extend(zip(*genome_data_addl)[1])
genome_names.extend(zip(*genome_data_addl)[3])

genome_dict = dict(zip(genome_ids, genome_names)) # NC_id: genome name
#pprint.pprint(genome_names)

# Get counts of each genome from fastq
picklename = 'fastq_counts.p'
if not os.path.exists(picklename):
	print "Counting fastq reads..."
	fastq = open(args.fastq_filename,'r')
	genome_counts = collections.Counter()
	for line in fastq:
		if line.startswith('@'):
			nc_id = line.partition('.fna')[0][1:] # @NC_007297.fna:1:1:1:5#0/1 -> NC_007297
			try:
				genome_counts.update([genome_dict[nc_id]]) # counter is keyed by genome name, not ID
			except:
				pass
				print line
	cPickle.dump(genome_counts,open(picklename,'wb'))
else:
	genome_counts = cPickle.load(open(picklename,'rb'))
print len(genome_counts)

#pprint.pprint(sorted(genome_counts.keys()))

with open('fastq_counts.csv', 'wb') as csvfile:
	counts_csv = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
	for k,v in genome_counts.items():
		counts_csv.writerow([k,v]) # csv writer needs a list, not a string

