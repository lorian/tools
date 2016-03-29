'''
many target_ids to 1 strain:
target_id strain
target_1 gene_1
target_2 gene_1
target_3 gene_2
'''
from M_compare_metagenome_results_to_truth import *
import csv
import argparse
import os
import cPickle

parser = argparse.ArgumentParser(description='Create t2g file for sleuth, linking individual fasta entries to their formal strain names.')
parser.add_argument('filename', help='Output file of kallisto')
args = parser.parse_args()
filename = args.filename

'''
global tax_dict
pickle_len = 0
if os.path.exists('species_taxonomy.pickle'):
	print "Loading taxonomy dict..."
	tax_dict = cPickle.load(open('species_taxonomy.pickle','rb'))
	pickle_len = len(tax_dict.names.keys())
'''
raw_est = Dataset()
input_csv = csv.reader(open(filename,'r'), 'excel-tab')
input_data = [r for r in input_csv]
input_data = input_data[1:] #remove header row
raw_est.species = zip(*input_data)[0]
raw_est.counts = raw_est.null_list()
raw_est.abundance = raw_est.null_list()

print "Number of raw entries: {0}".format(len(raw_est.species))

raw_est.remove_matches('rna') # remove specific genes
raw_est.remove_matches('gene_')
raw_est.remove_matches('plasmid')

with open("t2g_strains.txt","w") as f:
	f.write("target_id\tstrain\n")
	for st in raw_est.species:
		f.write("{}\t{}\n".format(st,st.partition('|')[0].partition('_gi')[0].replace('_',' ')))

raw_est.clean_names()
est = collapse_duplicates(raw_est)

'''
est.species = lookup_tax_list(est.species) # from this point on, species are taxids

if pickle_len != len(tax_dict.names.keys()): # don't re-pickle if nothing new is added
	print "Saving taxonomy dict..."
	cPickle.dump(tax_dict,open('species_taxonomy.pickle','wb'))
'''
