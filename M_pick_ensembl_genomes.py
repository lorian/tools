# Pick representative strains from ensembl bacteria dump

import os
import collections

def count_sp(fastas):
	# pick out unique species in directory files
	unique_species = collections.Counter()
	for f in fastas:
		name = f.partition('GCA')[0].partition('gca')[0].partition('.')[0].strip('_').strip('.')
		name_parts = name.split('_')
		if len(name_parts) == 2: # just the species is left
			unique_species.update([name])
		elif 'endosymbiont' in name_parts: # hard to parse these
			unique_species.update([name])
		elif '_sp._' in name or '_sp_' in name or '_species_' in name:
			unique_species.update(['_'.join(name_parts[0:3])])
		else:
			unique_species.update(['_'.join(name_parts[0:2])])
	return unique_species

def main():
	fastas = [d for d in os.listdir('.') if (d.endswith('dna.genome.fa') or d.endswith('cdna.all.fa')]
	unique_species = count_sp(fastas)

	print unique_species.most_common(50)
	print len(unique_species)
	thin_species = [u for u,c in unique_species.most_common(50) if c>100]
	fastas_delete = []
	for sp in thin_species: #remove half of each species' strains
		skipnext = False
		for f in fastas:
			if sp+"_" in f or sp+"." in f: # handle case where one sp is a longer version of other sp name
				if skipnext:
					fastas_delete.append(f)
					os.system("mv {} ignore/".format(f))
					skipnext = False
				else:
					skipnext = True

	print count_sp(fastas_delete)

if __name__ == "__main__":
	main()
