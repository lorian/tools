'''
The results are tab delimited lists of Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes:
Bacillus_aryabhattai_gca_001043825.GCA_001043825.1.30.dna.genome.fa     /home/lorian/scratch/illumina_100species_trimmed.1.fq   1       1       0/1000
Bacillus_atrophaeus_1942.GCA_000165925.1.30.dna.genome.fa       /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.295981        0.000944763     1/1000
Bacillus_bombysepticus_str_wang.GCA_000831065.1.30.dna.genome.fa        /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.197292        2.3647e-28      8/1000
'''

import csv
import argparse
import pprint
import collections
import os

def count_sp(fastas):
	# pick out unique species, more or less
	unique_species = collections.Counter()
	for f in fastas:
		name = f.partition('GCA')[0].partition('gca')[0].partition('.')[0].strip('_').strip('.')
		name_parts = name.split('_')
		if len(name_parts) == 2: # just the species is left
			sp_name = name
		elif 'endosymbiont' in name_parts: # hard to parse these
			sp_name = name
		elif '_sp._' in name or '_sp_' in name or '_species_' in name:
			sp_name = '_'.join(name_parts[0:3])
		else:
			sp_name = '_'.join(name_parts[0:2])

		if not (f.strip('_').startswith(sp_name+"_") or f.strip('_').startswith(sp_name+".")): # catch weird shit
			print "{} -> {} -> {}??".format(f,name,sp_name)

		unique_species.update([sp_name])
	return unique_species


parser = argparse.ArgumentParser(description='Process mash results for kallisto index creation')
parser.add_argument('filename', help='mash output file')
parser.add_argument('top_strains', help="How many strains of each species to keep for the quantification step")
args = parser.parse_args()

with open(args.filename,'r') as mash_file:
	mash_csv = csv.reader(mash_file, delimiter='\t')
	mash_data = [r for r in mash_csv]
	mash_1 = dict()
	mash_5 = dict()
	truth = dict()
	for r in mash_data:
		matching_hashes = int(r[4].split('/')[0])
		if r[0][0].islower():
			truth[r[0]] = matching_hashes
		if matching_hashes > 0:
			#print "{}\t\t{}\t{}".format(r[0],r[4],r[3])
			mash_1[r[0]] = matching_hashes
		if matching_hashes > 4:
			mash_5[r[0]] = matching_hashes

#pprint.pprint(truth)
species = count_sp(mash_5.keys())
sp_map = {sp.lower(): [(st,mash_1[st]) for st in mash_5.keys() if sp.lower() in st.lower()] for sp in species}
final_st = []
print args.top_strains
for sp in sp_map.keys():
	# pick top N of each species
	sp_map[sp].sort(key=lambda x:x[1],reverse=True)
	#print sp_map[sp][:min(args.top_strains,len(sp_map[sp]))]
	final_st.extend(sp_map[sp][:min(int(args.top_strains),len(sp_map[sp]))])
	#final_st.extend(sp_map[sp][:min(10,len(sp_map[sp]))])

#pprint.pprint(sp_map)
#pprint.pprint(truth)
#pprint.pprint(zip(*final_st)[0])
for f in zip(*final_st)[0]:
	os.system("mv {} ../".format(f))

print "All hits: {}\t At least 5 hits: {}\t Species: {}\t Strains kept: {}".format(len(mash_1.keys()), len(mash_5.keys()), len(species), len(final_st))