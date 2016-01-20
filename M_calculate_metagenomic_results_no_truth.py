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
import urllib
import urllib2
import pylab
import argparse
import copy

numpy.set_printoptions(precision=4)

class TaxaDict():
	# Consists of a dict mapping names -> taxids, and a dict mapping taxids -> Taxonomy objects
	def __init__(self, names = dict(), taxids = dict()):
		self.names = names
		self.taxids = taxids

	def get_taxa(self,name_or_id):
		taxa = self.get_taxa_by_id(name_or_id)
		if not taxa:
			taxa = self.get_taxa_by_name(name_or_id)
		return taxa

	def get_taxa_by_name(self,name):
		# get Taxonomy object by name
		try:
			return self.taxids[self.names[name]]
		except:
			try:
				return self.taxids[self.names[name.lower()]]
			except:
				return ''

	def get_taxa_by_id(self,taxid):
		# get Taxonomy object by ID
		try:
			return self.taxids[taxid]
		except:
			return ''

	def get_name_by_id(self,taxid):
		try:
			return self.taxids[taxid].name
		except:
			return ''

	def add_taxa(self,name,official_name,taxid,lineage):
		tax_entry = Taxonomy(official_name,taxid,lineage)
		self.names[str(name)] = taxid
		self.taxids[taxid] = tax_entry
		return tax_entry

	def remove_taxa(self,taxid):
		del self.names[self.get_taxa_by_id(taxid).name]
		del self.taxids[taxid]

# global taxa dictionary
tax_dict = TaxaDict()

class Taxonomy():

	def __init__(self, name = "", taxid = 0, lineage = dict()):
		self.name = name
		self.taxid = taxid
		self.lineage = lineage

	def __str__(self):
		return "{} | {} | {}".format(self.name,self.taxid,self.lineage)

	def filter_tier(self,rank):
		# return only if at or below a given taxonomic rank
		ranks = ['strain','species','genus','family','order','class','phylum']
		valid_ranks = set(ranks[0:ranks.index(rank)+1])
		if valid_ranks & set(self.lineage.keys()):
			return self.taxid
		else:
			return 0

	def filter_tier_exact(self,rank):
		# return only if AT a given taxonomic rank
		if rank in self.lineage.keys():
			return self.taxid
		else:
			return 0

	def return_tier(self,rank):
		# return the taxid for the specified rank
		try:
			return self.lineage[rank]
		except:
			return self.taxid

	def has_lineage(self,taxid):
		# does this Taxonomy object have this taxid somewhere in its lineage?
		return taxid in self.lineage.values()

class Dataset():

	def __init__(self, species=[], abundance=[], counts=[], size=[]):
		self.species = list(species)
		self.abundance = list(abundance)
		self.counts = list(counts)
		self.size = list(size) # can be empty

	def add_record(self, new_species, new_abundance, new_counts, new_size=None):
		self.species.extend([new_species])
		self.abundance.append(new_abundance)
		self.counts.append(new_counts)
		self.size.append(new_size)

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
		try:
			return self.counts[self.species.index(species)]
		except:
			return 0

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

		# Force clean_names to be lowercase and replace spaces with underscores
		self.species = [s.lower().strip().replace(" ","_") for s in clean_names]

	def sort_by_name(self):
		# Alphabetical by species
		self.set_by_array(sorted(self.get_array(),key=lambda x:x[0]))

	def set_threshold(self, threshold=0.001):
		# Set very low estimated abundances to 0
		self.set_by_array([r if r[1] > threshold else (r[0],0,r[2]) for r in self.get_array()])
		print "Number of non-zero abundances: {0}".format(len(self.species) - self.abundance.count(0))

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

def collapse_strains(strains,rank):
	""" Group strains together by taxa """
	just_rank = uptier_taxa(strains.species,rank)

	rank_combo = [x for x in zip(just_rank,strains.abundance,strains.counts,strains.size)]
	rank_dict = {}
	for s,a,c,z in rank_combo:
		rank_dict.setdefault(s,[0,0,0]) # stores ab,counts,size in that order
		rank_dict[s][0] = rank_dict[s][0] + a
		rank_dict[s][1] = rank_dict[s][1] + c
		rank_dict[s][2] = rank_dict[s][2] + z

	a,c,z = zip(*rank_dict.values())
	j_rank = Dataset(rank_dict.keys(),a,c,z)

	return j_rank

def collapse_duplicates(raw_data):
	# Create dictionary of lists of duplicates
	dup_data = raw_data.get_array()
	set_sp = {}
	set_ab = {}
	set_co = {}
	set_sz = {}
	set_plasmids = {}
	for sp,ab,co in dup_data:
		name = sp.partition('_gi|')[0] #the prepended strain name
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

	for k,v in set_sp.items():
		if len(v) == 1: # just add record directly if it has no duplicates
			undupe.add_record(k,set_ab[k][0],set_co[k][0]+set_plasmids[k],set_sz[k][0])
		else: # sum counts and average abundances
			undupe.add_record(k,math.fsum(set_ab[k])/len(v),math.fsum(set_co[k])+set_plasmids[k],math.fsum(set_sz[k]))

	print "Number of entries after combining duplicates: {0}".format(len(undupe.species))

	return undupe

def lookup_tax(original_name):
	global tax_dict

	known_names = dict([('bacterium ellin514 strain ellin514','Pedosphaera parvula Ellin514'),
						('bacteroidetes sp. f0058','Bacteroidetes oral taxon 274 str. F0058'),
						('baumannia cicadellinicola str. hc','Baumannia cicadellinicola str. Hc (Homalodisca coagulata)'),
						('bifidobacterium longum infantis atcc 55813','Bifidobacterium longum subsp. longum ATCC 55813'),
						('borrelia afzelii pko clone','Borrelia afzelii PKo'),
						('brucella abortus bv. 1 str. 9-941 chromosome','Brucella abortus bv. 1 str. 9-941'),
						('brucella abortus bv. 3 str.','Brucella abortus bv. 3'),
						('buchnera aphidicola str. bp','Buchnera aphidicola str. Bp (Baizongia pistaciae)'),
						('buchnera aphidicola str. sg','Buchnera aphidicola str. Sg (Schizaphis graminum)'),
						('buchnera aphidicola str. cc','Buchnera aphidicola BCc'),
						('buchnera aphidicola str. lsr1','Buchnera aphidicola str. LSR1 (Acyrthosiphon pisum)'),
						('candidatus hamiltonella defensa 5at','Candidatus Hamiltonella defensa 5AT (Acyrthosiphon pisum)'),
						('candidatus pelagibacter sp. htcc7211 1105874033148','Candidatus Pelagibacter sp. HTCC7211'),
						('candidatus ruthia magnifica str. cm','Candidatus Ruthia magnifica str. Cm (Calyptogena magnifica)'),
						('candidatus sulcia muelleri str. hc','Candidatus_Sulcia_muelleri_str._Hc_(Homalodisca_coagulata)'),
						('capnocytophaga sputigena atcc 33612 strain capno','Capnocytophaga sputigena ATCC 33612'),
						('clostridiales bacterium 1 7 47faa strain','Clostridiales bacterium 1_7_47FAA'),
						('clostridiales genomosp. bvab3 upii9-5','Mageeibacillus indolicus UPII9-5'),
						('clostridium difficile complete genome strain cf5','clostridium difficile cf5'),
						('clostridium difficile complete genome strain m120','clostridium difficile m120'),
						('enterococcus faecalis tx 0109','enterococcus faecalis tx0109'),
						('enterococcus faecalis tx 0411','enterococcus faecalis tx0411'),
						('enterococcus faecalis tx 0855','enterococcus faecalis tx0855'),
						('enterococcus faecalis tx 0860','enterococcus faecalis tx0860'),
						('enterococcus faecalis tx 2134','enterococcus faecalis tx2134'),
						('enterococcus faecalis tx 4248','enterococcus faecalis tx4248'),
						('escherichia coli o150:h5 se15','Escherichia coli SE15'),
						('lactobacillus brevis gravesensis atcc 27305','Lactobacillus brevis subsp. gravesensis ATCC 27305'),
						('lactobacillus reuteri sd2112 atcc 55730','Lactobacillus reuteri SD2112'),
						('lactobacillus rhamnosus gg atcc 53103','Lactobacillus rhamnosus GG'),
						('methanocaldococcus infernus me c','Methanocaldococcus infernus ME'),
						('nostoc azollae','\'Nostoc azollae\' 0708'),
						('providencia alcalifaciens ban1 integrating conjugative element icepalban1','Providencia alcalifaciens Ban1'),
						('salmonella enterica subsp. arizonae serovar 62:z4z23:--','Salmonella enterica subsp. arizonae serovar 62:z4,z23:-'),
						('salmonella enterica subsp. enterica serovar 4','Salmonella enterica subsp. enterica serovar 4,[5],12:i:- str. CVM23701'),
						('selenomonas sp. 67h29bp f0410','Selenomonas sp. oral taxon 149 str. 67H29BP'),
						('staphylococcus aureus aureus atcc baa-39','Staphylococcus aureus subsp. aureus ATCC BAA-39'),
						('staphylococcus aureus aureus mn8','Staphylococcus aureus subsp. aureus MN8'),
						('staphylococcus aureus aureus tch130/st-72','Staphylococcus aureus subsp. aureus TCH130'),
						('staphylococcus aureus aureus tch60','Staphylococcus aureus subsp. aureus TCH60'),
						('staphylococcus aureus aureus tch70','Staphylococcus aureus subsp. aureus TCH70'),
						('staphylococcus aureus aureus usa300 tch959','Staphylococcus aureus subsp. aureus USA300_TCH959'),
						('streptococcus gallolyticus tx20005','Streptococcus gallolyticus subsp. gallolyticus TX20005'),
						('streptococcus sp. 73h25ap f0408','Streptococcus sp. oral taxon 071 str. 73H25AP'),
						('streptomyces coelicolor a3','Streptomyces coelicolor A3(2)'),
						('wolbachia endosymbiont of drosophila willistoni tsc','Wolbachia endosymbiont of Drosophila willistoni TSC#14030-0811.24'),
						('xanthomonas campestris pv. vasculorum ncppb702','Xanthomonas vasicola pv. vasculorum NCPPB 702'),
						('bordetella bronchiseptica strain rb50','Bordetella bronchiseptica RB50'),
						('chlamydia trachomatis l2buch-1proctitis','Chlamydia trachomatis L2b/UCH-1/proctitis'),
						('ignicoccus hospitalis kin4i','Ignicoccus hospitalis strain KIN4/I'),
						('lawsonia intracellularis phemn1-00','Lawsonia intracellularis PHE/MN1-00'),
						('mycobacterium bovis af212297','Mycobacterium bovis AF2122/97'),
						('escherichia coli str. k-12 substr. mg1655star','Escherichia coli str. K-12 substr. MG1655'),
						('clostridium botulinum e1 str. \'bont e','Clostridium botulinum E1 str. \'BoNT E Beluga\''),
						('bacteria',2), # taxa with non-bacterial entries in NCBI database
						('bacillus',1386),
						('thermosipho',2420),
						('yersinia',629),
						])

	if original_name.lower().replace('_',' ').strip() in known_names.keys(): # problematic names
		name = known_names[original_name.lower().replace('_',' ').strip()]
	else:
		name = original_name

	tax_entry = tax_dict.get_taxa_by_name(original_name) # skip the lookup if taxa was already found
	if tax_entry:
		return tax_entry
	tax_entry = tax_dict.get_taxa_by_id(name) # skip the lookup if taxa was already found
	if tax_entry:
		return tax_entry

	try:
		if type(name) == int or name.isdigit(): #if we've already passed a taxid instead of a text name
			taxid = name
		else:
			# Clean name for URL use; replace marks that cause problems for NCBI lookup with spaces
			url_name = urllib.quote_plus(name.translate(string.maketrans("()[]:","     ").replace('_',' ').strip()))

			# Look up taxonomy ID
			handle = Entrez.esearch(db="taxonomy", term=url_name)
			records = Entrez.read(handle)
			taxid = records['IdList']
			if not taxid: # lookup manually instead of through biopython

				page_taxa = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={}'.format(url_name))
				for line in page_taxa:
					line = string.replace(line,"\t","")
					if line.startswith("<Id>"):
						taxid = line[4:-6]

				if not taxid:
					print "\t Unable to find {}".format(name)
					print url_name
					return ""

		if str(taxid) == '178505':
			print original_name
			print url_name


		# Get taxonomy for id
		handle = Entrez.efetch(db="taxonomy", id=taxid, mode="text", rettype="xml")
		taxon = Entrez.read(handle)[0] # only grab the first

		taxid = taxon["TaxId"]
		official_name = taxon['ScientificName']
		lineage = dict()
		try:
			for t in taxon["LineageEx"]:
				lineage[t['Rank']] = t['TaxId']
		except:
			pass # no parent taxonomy

		if taxon['Rank'] == 'no rank':
			if 'species' in lineage.keys():
				lineage['strain'] = taxid
			else:
				lineage[str(taxid)] = taxid
		else:
			lineage[taxon['Rank']] = taxid

		if type(name) == int or name.isdigit():
			tax_entry = tax_dict.add_taxa(official_name,official_name,taxid,lineage)
			tax_entry = tax_dict.add_taxa(original_name,official_name,taxid,lineage)
			print "{} -> {}".format(original_name,official_name)
		else:
			tax_entry = tax_dict.add_taxa(original_name,official_name,taxid,lineage)
			print "{} -> {}".format(original_name,official_name)

		if '2759' in lineage.values(): # catch eukaryotes
			print "\t\t\tWarning: {} appears to be a eukaryote!".format(original_name)
			print tax_entry

	except Exception as e: # because when a long string of name lookups errors in the middle, it hurts
		print e
		return ""

	return tax_entry

def lookup_tax_list(species_list):
	species_tax = []
	for name in species_list:
		tax_entry = lookup_tax(name)
		if tax_entry:
			species_tax.append(tax_entry.taxid)

	return species_tax

def process_input(filename):
	""" Pull species names, abundance, and counts out of input file """

	suffix = filename.rpartition('.')[2]
	input_file = open(filename,'r')

	raw_est = Dataset()

	input_csv = csv.reader(input_file, 'excel-tab')
	input_data = [r for r in input_csv]
	input_data = input_data[1:] #remove header row
	raw_est.species = zip(*input_data)[0]
	raw_est.counts = [float(i) for i in zip(*input_data)[3]]
	raw_est.abundance = [float(i) for i in zip(*input_data)[4]]
	raw_est.size = raw_est.null_list()

	print "Number of raw entries: {0}".format(len(raw_est.species))

	raw_est.clean_names()
	raw_est.remove_matches('rna') # remove specific genes
	raw_est.remove_matches('gene_')
	raw_est.remove_matches('plasmid')

	est = collapse_duplicates(raw_est)
	est.convert_to_percentage()
	est.sort_by_name()
	est.set_threshold()
	est.species = lookup_tax_list(est.species) # from this point on, species are taxids

	return est

def filter_by_taxa(taxids,rank):
	filtered_names = set()
	for taxid in taxids:
		id_taxa = tax_dict.get_taxa_by_id(taxid)
		filtered_names.add(id_taxa.filter_tier(rank))

	filtered_names.discard(0) # 0 is returned if a name doesn't pass the rank filter

	# Discard specific entries that don't get caught automatically
	if rank == 'strain':
		filtered_names.discard('83333') # Escherichia coli K-12

	return filtered_names

def climb_tree(taxid,dataset):
	# Given a taxid and a Dataset, calculate the total counts that would end up in that taxid
	count_sum = 0.0
	for sp in dataset.species:
		if tax_dict.get_taxa_by_id(sp).has_lineage(taxid):
			count_sum += dataset.lookup_count(sp)
			#print tax_dict.get_taxa_by_id(sp)
	return count_sum

def climb_tree_verbose(taxid,dataset):
	# Given a taxid and a Dataset, calculate the total counts that would end up in that taxid
	count_sum = 0.0
	for sp in dataset.species:
		if tax_dict.get_taxa_by_id(sp).has_lineage(taxid):
			count_sum += dataset.lookup_count(sp)
			print tax_dict.get_taxa_by_id(sp)
	return count_sum

def graph_est(est, expname, tier, show_graphs, save_graphs):
	# These imports are here so script can run on server w/out graphics
	import lanthplot
	import matplotlib
	import seaborn

	# Drop any entry that is above the targeted taxa level
	filtered_names = list(filter_by_taxa(est.species,tier))

	min_abundances = 0.01
	filtered_est = Dataset()
	for sp in filtered_names:
		if est.lookup_abundance(sp) > min_abundances: # filter out results with less than min abundances
			filtered_est.add_record(sp,est.lookup_abundance(sp),est.lookup_abundance(sp),est.lookup_size(sp))

	present_species = [tax_dict.get_name_by_id(s) for s in filtered_est.species]
	present_est = list(filtered_est.abundance)

	# graph est abundances
	xmax = len(present_species)
	present_sp = [x.replace('_',' ') for x in present_species]
	x = numpy.array(range(0,xmax))

	# sort based on est abundances
	all_filter = zip(present_sp,present_est)
	all_filter.sort( key=lambda x: x[1],reverse=True )
	fil_sp,fil_est = zip(*all_filter)

	lanthplot.plot_setup_pre(
		"Graph of all above-average estimated abundances in {} at {}-level"
		.format(expname,tier), xlabels = fil_sp, xticks = range(0,xmax),
		xrotation = -90, yaxislabel = 'Abundance')

	lanthplot.plot(x, fil_est, color='red', label="Estimated", plot_type = 'scatter')
	matplotlib.pyplot.gca().set_ylim(bottom=0.)
	if save_graphs:
		lanthplot.plot_setup_post(save_file = expname +'_'+ tier +'_estabundances.png', show=show_graphs)
	else:
		lanthplot.plot_setup_post(legend=False)

	return

def main(argv=sys.argv):
	"""
	Command line usage: python compare_metagenomic_results_to_truth.py
						[filename] [program (default kallisto)] <-g show graphs?>
						<-s save graphs?> <--taxa level (defaults to all)>
	"""

	parser = argparse.ArgumentParser(description='Calculates and displays metagenomic abundance, as estimated by kallisto')
	parser.add_argument('filename', help='Output file of metagenomic analysis tool')

	parser.add_argument('-g', '--show-graphs', action='store_true', help='Display graphs of calculated errors')
	parser.add_argument('-s','--save-graphs', action='store_true', help='Save graphs of calculated errors to file')
	parser.add_argument('--taxa', default='all', help='Desired taxa level of analysis. Accepts one of: strain, species, genus, phylum. Defaults to all levels.')
	args = parser.parse_args()

	filename = args.filename
	exp_name = filename.rpartition('.')[0] # will be used for graph-naming purposes
	show_graphs = args.show_graphs
	save_graphs = args.save_graphs

	print "Running analysis on {}\n".format(filename)

	global tax_dict
	pickle_len = 0
	if os.path.exists('species_taxonomy.pickle'): # stored taxonomy information from previous runs
		print "Loading taxonomy dict..."
		tax_dict = cPickle.load(open('species_taxonomy.pickle','rb'))
		pickle_len = len(tax_dict.names.keys())

	estimated = process_input(filename)
	est_j_species = collapse_strains(estimated,'species')
	est_j_genus = collapse_strains(estimated,'genus')
	est_j_phylum = collapse_strains(estimated,'phylum')

	if pickle_len != len(tax_dict.names.keys()): # don't re-pickle if nothing new is added
		print "Saving taxonomy dict..."
		cPickle.dump(tax_dict,open('species_taxonomy.pickle','wb'))

	dataset_pairs = [('strain',estimated),('species',est_j_species),('genus',est_j_genus),('phylum',est_j_phylum)]

	if args.taxa == 'strain':
		dataset_pairs = [('strain',estimated)]
	elif args.taxa == 'species':
		dataset_pairs = [('species',est_j_species)]
	elif args.taxa == 'genus':
		dataset_pairs = [('genus',est_j_genus)]
	elif args.taxa == 'phylum':
		dataset_pairs = [('phylum',est_j_phylum)]

	for label,est in dataset_pairs:
		graph_est(est, exp_name, label, show_graphs, save_graphs)

if __name__ == "__main__":
	main()
