# Process large metagenomic reference dataset:
#	Make names fully readable by replacing spaces with _ (done)
#	Discard virus genomes (done)
#	Combine all chromosomal entries for a single species using 10 N's (done)
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

# Open files
fa_genomes = open('genome_processed.fa', 'w')
fa_genomes.seek(0) #overwrite if it exists
fa_plasmids = open('plasmids_processed.fa', 'w')
fa_plasmids.seek(0) #overwrite if it exists

last_species = ""
last_type = ""

for line in mfa:
	if line[0] == '>': # fasta name
		line = line.replace (" ", "_")
		if line.find('VIRL') != -1: #ignore viruses
			last_type = 'virus'
		elif line.find('plasmid') != -1: #just copy plasmids
			last_type = 'plasmid'
			fa_plasmids.write(line)
		else: #attempt to "assemble" genomes
			last_type = 'genome'
			species,x,y = line.partition('|')
			if last_species == species:
				fa_genomes.write('NNNNNNNNNN') # potential gap
			else:
				fa_genomes.write(line)
				last_species = species
	else: #actual fasta content
		if last_type == 'genome':
			fa_genomes.write(line)
		elif last_type == 'plasmid':
			fa_plasmids.write(line)

# close files
fa_genomes.truncate()
fa_genomes.close()
fa_plasmids.truncate()
fa_plasmids.close()
