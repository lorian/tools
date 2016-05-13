"""
Compares count-based output of metagenomic analysis tools with known ground truth of dataset
Lorian Schaeffer, Pachter Lab, UC Berkeley
"""

import sys
import csv
import numpy
import string
import math
import itertools
import collections
import cPickle
import os
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'
import urllib
import urllib2
import argparse
import copy
import pprint

numpy.set_printoptions(precision=4)


class Dataset():

	def __init__(self, species=[], abundance=[], counts=[], size=[]):
		self.species = list(species)
		self.abundance = list(abundance)
		self.counts = list(counts)
		self.size = list(size) # can be empty
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts)) # defaults to 0

	def add_record(self, new_species, new_abundance, new_counts, new_size=None):
		self.species.extend([new_species])
		self.abundance.append(new_abundance)
		self.counts.append(new_counts)
		self.size.append(new_size)
		self.counts_by_sp[new_species] = new_counts

	def set_by_array(self, array):
		# Warning: this overwrites all existing data in this dataset
		try:
			self.species, self.abundance, self.counts, self.size = zip(*array)
		except:
			self.species, self.abundance, self.counts = zip(*array)
		self.species = list(self.species)
		self.abundance = [float(x) for x in self.abundance]
		self.counts = [float(x) for x in self.counts] # some programs like eXpress return counts as floats
		self.size = [int(x) for x in self.size]
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts))

	def remake_index(self):
		# when the species listing has updated and the dictionary needs to match
		self.counts_by_sp = collections.defaultdict(int,zip(self.species,self.counts))

	def print_record(self, species):
		print "{0}:\n\t{1} ({2} reads)".format(species,self.lookup_abundance(species),self.lookup_count(species))

	def get_array(self):
		return zip(self.species,self.abundance,self.counts)

	def lookup_abundance(self, species):
		try:
			return self.abundance[self.species.index(species)]
		except:
			for sp in species.split('?'):
				try:
					return self.abundance[self.species.index(sp)]
				except:
					pass
		return 0

	def lookup_count(self, species):
		return self.counts_by_sp[species]

	def lookup_size(self, species):
		try:
			return self.size[self.species.index(species)]
		except:
			return 0

	def null_list(self):
		# For use when one of the variables needs to be filled in
		return [0]*len(self.species)

	def clean_names(self):
		# If BACT_ abbreviations are found, see if they match large genome database
		'''
		try:
			min( i for i, sp in enumerate(self.species) if 'BACT_' in sp )
		except:
			clean_names = self.species
		else:
			with open('fullgenomenames.txt', 'r') as biggenomefile:
				big_genome = dict(csv.reader(biggenomefile))
			clean_names = [big_genome[s.partition('|')[0]]+"_"+s.partition('|')[2]
							if s.partition('|')[0] in big_genome.keys()
							else s for s in self.species]
		'''
		# Force clean_names to be lowercase and replace spaces with underscores
		self.species = [s.lower().strip().replace(" ","_") for s in self.species]

	def sort_by_name(self):
		# Alphabetical by species
		self.set_by_array(sorted(self.get_array(),key=lambda x:x[0]))

	def set_threshold(self, threshold=1):
		# Drop entries with very low estimated counts
		self.set_by_array([r for r in self.get_array() if r[2] >= threshold])
		print "Number of non-zero abundances: {0}".format(len(self.species))

	def convert_to_percentage(self):
		# Normalize TPM/FPKM abundances
		total_ab = math.fsum(self.abundance)
		if total_ab == 0: # give up if there are no abundances
			return

		unknown = self.lookup_abundance('unknown')
		self.abundance = [100*v/(total_ab - unknown) for v in self.abundance]
		if unknown != 0: # prevent unknowns from being normalized
			self.abundance[self.species.index("unknown")] = unknown

		assert round(math.fsum(self.abundance)) == 100 + round(unknown)

	def remove_matches(self, target):
		self.set_by_array([r for r in self.get_array() if (r[0].find(target) == -1)])
		print "Number of entries after removing {0}: {1}".format(target,len(self.species))

	def delete_record(self,target):
		try:
			index = self.species.index(target)
		except:
			print "{} not found to delete".format(target)
			return
		del self.species[index]
		del self.counts[index]
		del self.abundance[index]
		del self.counts_by_sp[index]
		try:
			del self.size[index]
		except:
			return


def print_names(strains):
	print [tax_dict.get_name_by_id(s) for s in strains]

def uptier_taxa(strains,rank):
	uptier_names = filter_by_taxa(strains,rank) # gives entries of rank or below

	uptier_ids = []
	for s in strains:
		if s in uptier_names:
			rank_id = tax_dict.get_taxa_by_id(s).return_tier(rank)
			uptier_ids.append(rank_id)
		else:
			uptier_ids.append(s)

	lookup_tax_list(uptier_ids)
	return uptier_ids


def collapse_duplicates(raw_data):
	# Create dictionary of lists of duplicates
	dup_data = raw_data.get_array()
	set_sp = {}
	set_ab = {}
	set_co = {}
	set_sz = {}
	set_plasmids = {}
	for sp,ab,co in dup_data:
		name = sp.partition('_gi|')[0].partition('|')[0] #the prepended strain name
		set_plasmids.setdefault(name,0) # so there's always a plasmid count value for any given name key
		if 'plasmid' in sp and not (name == 'ralstonia_eutropha_h16' or name == 'cupriavidus_necator_h16'):
			set_plasmids[name] += co
		else:
			set_sp.setdefault(name,[]).append(sp)
			set_ab.setdefault(name,[]).append(ab)
			set_co.setdefault(name,[]).append(co)
			set_sz.setdefault(name,[]).append(co)

	assert(set_ab.keys() == set_co.keys() == set_sp.keys())

	# New, clean dataset for data without duplicates
	undupe = Dataset()

	# Note: we include plasmids in the count total solely because i100 was simulated to include 1x plasmid coverage.
	for k,v in set_sp.items():
		if len(v) == 1: # just add record directly if it has no duplicates
			undupe.add_record(k,set_ab[k][0],set_co[k][0]+set_plasmids[k],set_sz[k][0])
		else: # sum counts and average abundances
			undupe.add_record(k,math.fsum(set_ab[k])/len(v),math.fsum(set_co[k])+set_plasmids[k],math.fsum(set_sz[k]))

	print "Number of entries after combining duplicates: {0}".format(len(undupe.species))

	return undupe

def fix_transcript_names(species):
	cleanup = string.maketrans('-()+/','_____')
	transcript_names = [(sp.partition("|")[0].translate(cleanup,'.[]') +"_"+
		"..".join([d for d in sp.split(':') if d.isdigit() and len(d)>1])) for sp in species]

	return transcript_names


def process_input(filename,program,truth,fragmented=False):
	""" Pull species names, abundance, and counts out of input file """

	suffix = filename.rpartition('.')[2]

	if os.path.exists(filename +'.p'): # kraken output is slow to process, so look for saved processed version
		est = cPickle.load(open(filename +'.p','rb'))
		est.set_threshold()
		est.species = lookup_tax_list(est.species) # from this point on, species are taxids
		return est
		#input_table = cPickle.load(open(filename +'.p','rb'))
	else:
		input_file = open(filename,'r')

	raw_est = Dataset()
	raw_est.size = truth.size
	if program == 'express':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[1]
		raw_est.clean_names()
		raw_est.counts = [float(i) for i in zip(*input_data)[6]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[10]]
		raw_est.size = raw_est.null_list()

	elif program == 'clark':
		input_csv = csv.reader(input_file, 'excel')
		input_data = [r for r in input_csv]
		input_data = input_data[1:-1] #remove header row and unknown line at end
		raw_est.species = zip(*input_data)[0]
		raw_est.clean_names()
		raw_est.counts = [float(i) for i in zip(*input_data)[2]]
		raw_est.abundance = [float(i) if i != '-' else 0.0 for i in zip(*input_data)[4]]

	elif program == 'kraken':
		if not os.path.exists(filename +'.p'):
			input_csv = csv.reader(input_file, 'excel-tab')
			input_counts = collections.Counter()
			for r in input_csv:
				input_counts.update([r[1].rpartition(';')[2]])
			input_table = [(k.rpartition(';')[2],v) for k,v in input_counts.iteritems()]
			cPickle.dump(input_table,open(filename +'.p','wb')) # save processed form
		raw_est.species = zip(*input_table)[0]
		raw_est.counts = [float(i) for i in zip(*input_table)[1]]
		raw_est.abundance = [0]*len(raw_est.species)

	elif program == 'kallisto':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[0]
		raw_est.counts = [float(i) for i in zip(*input_data)[3]]
		raw_est.abundance = [float(i) for i in zip(*input_data)[4]]
		raw_est.size = raw_est.null_list()

	elif program == 'gasic':
		input_csv = csv.reader(input_file, 'excel-tab')
		input_data = [r for r in input_csv]
		input_data = input_data[1:] #remove header row
		raw_est.species = zip(*input_data)[0]
		raw_est.counts = [float(i) for i in zip(*input_data)[2]]
		raw_est.abundance = [0]*len(raw_est.species)

	else:
		print "File is not supported input type"

	print "Number of raw entries: {0}".format(len(raw_est.species))

	raw_est.clean_names()
	if not transcripts:
		raw_est.remove_matches('rna') # remove specific genes
		raw_est.remove_matches('gene_')
		est = collapse_duplicates(raw_est)
	else:
		raw_est.remove_matches('plasmid')
		raw_est.species = fix_transcript_names(raw_est.species)
		est = raw_est

	est.convert_to_percentage()
	est.sort_by_name()
	est.set_threshold()

	cPickle.dump(est,open(filename +'.p','wb'))

	est.remake_index()
	return est

def main(argv=sys.argv):
	"""
	Command line usage: python compare_metagenomic_results_to_truth.py
						[filename] [program (default kallisto)] <-g show graphs?>
						<-s save graphs?> <--taxa level (defaults to all)>
	"""

	parser = argparse.ArgumentParser(description='Compare output of metagenomic analysis tools with ground truth of dataset')
	parser.add_argument('filename', help='Output file of metagenomic analysis tool')
	parser.add_argument('program', nargs='?', default='kallisto', help='Source program that created output. Valid options are: kallisto, kraken, clark, gasic, express. Defaults to kallisto.')

	parser.add_argument('-g', '--show-graphs', action='store_true', help='Display graphs of calculated errors')
	parser.add_argument('-s','--save-graphs', action='store_true', help='Save graphs of calculated errors to file')
	parser.add_argument('--taxa', default='all', help='Desired taxa level of analysis. Accepts one of: strain, species, genus, phylum. Defaults to all levels.')
	parser.add_argument('--dataset', default='i100', help='Dataset truth to be compared to. Accepts: i100, simmt_have, simmt_all, simmt_transcripts, no_truth. Defaults to i100.')

	args = parser.parse_args()

	filename = args.filename
	exp_name = filename.rpartition('.')[0] # will be used for graph-naming purposes
	program = args.program
	show_graphs = args.show_graphs
	save_graphs = args.save_graphs
	dataset = args.dataset

	global transcripts
	transcripts = False
	if dataset == 'simmt_transcripts':
		transcripts = True

	print "Running comparison on {}\n".format(filename)

	global tax_dict
	pickle_len = 0
	if os.path.exists('species_taxonomy.pickle') and not transcripts: # stored taxonomy information from previous runs
		print "Loading taxonomy dict..."
		tax_dict = cPickle.load(open('species_taxonomy.pickle','rb'))
		pickle_len = len(tax_dict.names.keys())

	estimated = process_input(filename,program,truth)

	if pickle_len != len(tax_dict.names.keys()): # don't re-pickle if nothing new is added
		print "Saving taxonomy dict..."
		cPickle.dump(tax_dict,open('species_taxonomy.pickle','wb'))

if __name__ == "__main__":
	main()
