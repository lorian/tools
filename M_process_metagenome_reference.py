# Process large metagenomic reference dataset:
#	Make names fully readable by replacing spaces with _ (done)
#	Discard virus genomes (done)
#	Combine all chromosomal entries for a single species using 10 N's
#	Assemble plasmids??
#	Split into files of less than 3.6 billion characters for bowtie2 index

import sys

filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if filename == "":
		filename = arg #avoid extra space at beginning
	else:
		filename = filename + " " + arg

mfa = open(filename,'r')
text = ""
plasmids = ""
genomes = dict([])
neither = ""
last_species = ""

for line in mfa:
	if line.startswith('>'):
		text = text + line.replace (" ", "_")
		if line.find('VIRL') != -1: #ignore viruses
			pass
		elif line.find('plasmid') != -1: #just copy plasmids
			plasmids = plasmids + line
		elif line.find('genome') != -1: #attempt to "assemble" genomes
			species,x,y = line.partition('|')
			if last_species == species:
				genomes[species] = genomes[species] + 'NNNNNNNNNN' + line
			else:
				genomes[species] = line
				last_species = species
		else:
			neither = neither + line
#	else:
#		text = text + line

# copy all fasta names to file
fa_genomes = open('genome_names.txt', 'w')
fa_genomes.seek(0) #overwrite if it exists
for g in genomes:
	fa_genomes.write(g + '\n' + genomes[g] + '\n')
fa_genomes.truncate()
fa_genomes.close()

fa_plasmids = open('plasmid_names.txt', 'w')
fa_plasmids.seek(0) #overwrite if it exists
fa_plasmids.write(plasmids)
fa_plasmids.truncate()
fa_plasmids.close()

fa_other = open('other_names.txt', 'w')
fa_other.seek(0) #overwrite if it exists
fa_other.write(neither)
fa_other.truncate()
fa_other.close()

#fa = open(filename[:filename.rfind('.')] + '_processed.fa', 'w')
#fa.seek(0) #overwrite if it exists
#fa.write(text)
#fa.truncate()
#fa.close()
