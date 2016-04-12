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

def main():
	fastas = [d for d in os.listdir('.') if (d.endswith('dna.genome.fa') or d.endswith('cdna.all.fa'))]
	unique_species = count_sp(fastas)

	fastas_delete = set()
	print "{} species".format(len(unique_species))
	'''
	thin_species = [u for u,c in unique_species.most_common(50) if c>100]

	for sp in thin_species: #remove half of each species' strains
		skipnext = False
		for f in fastas:
			if sp+"_" in f or sp+"." in f: # avoid case where one sp is a longer version of other sp name
				if skipnext:
					fastas_delete.add(f)
					os.system("mv {} ignore/".format(f))
					skipnext = False
				else:
					skipnext = True
	'''
	dup_species = [u for u,c in unique_species.most_common() if c>1]
	for sp in dup_species: #only keep largest strain of each species
		fastas_size = []
		for f in fastas:
			if f.strip('_').startswith(sp+"_") or f.strip('_').startswith(sp+"."): # avoid case where one sp is a longer version of other sp name
				fastas_size.append((os.path.getsize(f),f))
				fastas_delete.add(f)
		sizes,names = zip(*fastas_size)
		big_sp = names[sizes.index(max(sizes))]
		fastas_delete.discard(big_sp)
		assert os.path.getsize(big_sp) == max(sizes)

	print "{} total fasta files".format(len(fastas))
	print "{} fastas moved to /ignore".format(len(fastas_delete))
	print "{} fastas remaining".format(len(fastas) - len(fastas_delete))
	for f in fastas_delete:
		os.system("mv {} ignore/".format(f))
	#print count_sp(fastas_delete)

if __name__ == "__main__":
	main()
