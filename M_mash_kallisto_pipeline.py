'''
The results are tab delimited lists of Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes:
Bacillus_aryabhattai_gca_001043825.GCA_001043825.1.30.dna.genome.fa     /home/lorian/scratch/illumina_100species_trimmed.1.fq   1       1       0/1000
Bacillus_atrophaeus_1942.GCA_000165925.1.30.dna.genome.fa       /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.295981        0.000944763     1/1000
Bacillus_bombysepticus_str_wang.GCA_000831065.1.30.dna.genome.fa        /home/lorian/scratch/illumina_100species_trimmed.1.fq   0.197292        2.3647e-28      8/1000
'''

# run from inside directory with ALL genomes; will move genomes one level up
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

def collapse_contigs(f):
	basename = f.partition('.dna.genome.fa')[0].partition('.mfa')[0]
	try:
		mfa = open(f,'r')
	except IOError:
		print "{} does not exist".format(f)
		return False

	text = ""

	firstline = True
	for line in mfa:
		if firstline:
			# replace header with name
			new_name = basename.partition('.1.30')[0].partition('.GCA')[0].replace(" ", "_")
			text = '>' + new_name + "|\n"
			firstline = False
		elif line.startswith('>'):
			text += 'NNNNNNNNNN' #indicate possible gaps between chr, plasmids, shotgun pieces, etc
		elif line!= "\n": #skip empty lines within fasta
			text += line

	text = text + "\n" #add newline to end of file for eventual cat
	fa = open(basename + '.cat.fa', 'w')
	fa.seek(0)
	fa.write(text)
	fa.truncate()
	fa.close()
	return basename + '.cat.fa'

parser = argparse.ArgumentParser(description='Process mash results for kallisto index creation')
parser.add_argument('filename', help='mash output file')
parser.add_argument('top_strains', help="How many strains of each species to keep for the quantification step")
parser.add_argument('directory', default="../", help="Directory to put files for kallisto index creation")
parser.add_argument('dry-run', default=False, help="If True, outputs files that would be moved, but does not move any files.")
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

for sp in sp_map.keys():
	# pick top N of each species
	sp_map[sp].sort(key=lambda x:x[1],reverse=True)
	#print sp_map[sp][:min(args.top_strains,len(sp_map[sp]))]
	final_st.extend(sp_map[sp][:min(int(args.top_strains),len(sp_map[sp]))])

final_names = list(set(zip(*final_st)[0]))
final_names.sort()
with open('mash_names.txt','w') as f:
	f.writelines(final_names)

if not args.dry_run: #don't create or move files

	for f in final_names:
		new_name = collapse_contigs(f)
		if new_name:
			os.system("mv {0} {1}".format(new_name,args.directory))
			if not os.path.exists(os.path.join('{0}','{1}').format(args.directory,new_name)):
				print new_name

	#compare lists
	moved_files = [f for f in os.listdir(args.directory) if f.endswith(".cat.fa")]
	for f in final_names:
		basename = f.partition('.dna.genome.fa')[0].partition('.mfa')[0]
		if not basename+'.cat.fa' in moved_files:
			print "Missing: {}".format(f)

	print "All hits: {}\t At least 5 hits: {}\t Species: {}\t Strains kept: {}\t Strains moved: {}\t".format(len(mash_1.keys()), len(mash_5.keys()), len(species), len(final_names), len(moved_files))

else:
	print "All hits: {}\t At least 5 hits: {}\t Species: {}\t Strains kept: {}\t".format(len(mash_1.keys()), len(mash_5.keys()), len(species), len(final_names))
