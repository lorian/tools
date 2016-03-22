# Add strain name to entries in refseq multifastas, and replace spaces with underscores

import os
import fileinput
import csv
from Bio import Entrez
Entrez.email = 'lanthala@berkeley.edu'
import urllib
import urllib2
import string
import cPickle

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

def lookup_tax(ncbi_ref,raw_name):
	global tax_dict

	tax_entry = tax_dict.get_taxa_by_name(ncbi_ref) # skip the lookup if taxa was already found
	if tax_entry:
		print tax_entry
		return tax_entry

	try:
		if type(ncbi_ref) == int or ncbi_ref.isdigit(): #if we've already passed a taxid instead of a text name
			taxid = ncbi_ref
		else:
			# Look up taxonomy ID
			page_taxa = urllib2.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=taxonomy&id={}'.format(ncbi_ref))
			nextline = False
			for line in page_taxa:
				line = line.strip()
				if line.startswith("<Link>"):
					nextline = True
				if nextline and line.startswith("<Id>"):
					taxid = line[4:-5]
					break

			if not taxid:
				print "\t Unable to find {}".format(ncbi_ref)
				print ncbi_ref
				return ""

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

		if 'species' in lineage.keys():
			if taxon['Rank'] == 'no rank':
				lineage['strain'] = taxid
			else:
				lineage[taxon['Rank']] = taxid
		else: # because "type strains" aren't marked as such
			lineage[taxon['Rank']] = taxid
			if ('strain' in raw_name or 'str' in raw_name) and 'str' not in official_name:
				if 'strain' in raw_name:
					splitname = raw_name.partition('strain_')
				elif 'str' in raw_name:
					splitname = raw_name.partition('strain_')
				print "\tERROR1: {}".format(official_name)
				official_name = (splitname[0] + splitname[1] + splitname[2].partition('_')[0]).strip('_')

		# Handle the case where the fasta name isn't the official name
		if official_name.replace(' ','_') not in raw_name:
			print "\tERROR2: {}".format(official_name)
			official_name = raw_name.partition('contig')[0].partition('Contig')[0].partition('chr')[0].partition('Chr')[0].strip('_')


		if type(ncbi_ref) == int or ncbi_ref.isdigit():
			tax_entry = tax_dict.add_taxa(official_name,official_name,taxid,lineage)
			tax_entry = tax_dict.add_taxa(ncbi_ref,official_name,taxid,lineage)
			print "{} -> {}".format(ncbi_ref,official_name)
		else:
			tax_entry = tax_dict.add_taxa(ncbi_ref,official_name,taxid,lineage)
			print "{} -> {}".format(ncbi_ref,official_name)

		if '2759' in lineage.values(): # catch eukaryotes
			print "\t\t\tWarning: {} appears to be a eukaryote!".format(ncbi_ref)
			print tax_entry

	except Exception as e: # because when a long string of name lookups errors in the middle, it hurts
		print e
		return ""

	return official_name

def main():
	global tax_dict
	unique_species = set()
	unique_strains = set()

	pickle_len = 0
	if os.path.exists('ncbi_ref_taxonomy.pickle'): # stored taxonomy information from previous runs
		print "Loading taxonomy dict..."
		tax_dict = cPickle.load(open('ncbi_ref_taxonomy.pickle','rb'))
		pickle_len = len(tax_dict.names.keys())
		print "Loading strain dict..."
		unique_strains = cPickle.load(open('ncbi_ref_strains.pickle','rb'))

	# Process fasta entries
	headers = open('bacteria_headers.txt','r')
	header_cols = csv.reader(headers,delimiter='|')
	for line in header_cols:
		if not line[4].lstrip('_').startswith(tuple(unique_strains)):
			#print line
			unique_strains.add(lookup_tax(line[3],line[4]).replace(' ','_'))

	if pickle_len != len(tax_dict.names.keys()): # don't re-pickle if nothing new is added
		print "Saving taxonomy dict..."
		cPickle.dump(tax_dict,open('ncbi_ref_taxonomy.pickle','wb'))
		print "Saving unique strains..."
		cPickle.dump(unique_strains,open('ncbi_ref_strains.pickle','wb'))

if __name__ == "__main__":
	main()

