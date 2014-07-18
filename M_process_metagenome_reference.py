# Process large metagenomic reference dataset:
#	Make names fully readable by replacing spaces with _ (done)
#	Discard virus genomes (done)
#	Combine all chromosomal entries for a single species using 10 N's (done)
#	Assemble plasmids?? (no)
#	Split into files of less than 3.6 billion characters for bowtie2 index (done)
#		note that the plasmid file never gets too big, so I don't bother checking if I should split it

#	check that all lines begin as expected: with > or with nucleotides

import sys
import os

split_size = 3500000000 # 3.5gb

#Note: File needs to be sorted!
filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if filename == "":
		filename = arg #avoid extra space at beginning
	else:
		filename = filename + " " + arg

mfa = open(filename,'r')

file_count_g = 0
file_count_p = 0
line_count = 0

# Open files
genomes_string = filename.rsplit('.', 1)[0] + "_genomes_{0}.mfa" #mfa file format is for compatibility with my mass-indexing script
plasmids_string = filename.rsplit('.', 1)[0] + "_plasmids_{0}.mfa"

fa_genomes = open(genomes_string.format(file_count_g), 'wt')
fa_genomes.seek(0) #overwrite if it exists
fa_plasmids = open(plasmids_string.format(file_count_p), 'wt')
fa_plasmids.seek(0) #overwrite if it exists

last_species = "."
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
				print "Species duplicate of {0} in {1}".format(species,line)
				fa_genomes.write('NNNNNNNNNN') # potential gap
			else:
				print "Species NOT duplicate of {0} in {1}".format(last_species,line)
				# Split file past size limit, but only between species entries
				fa_genomes.flush()
				filesize = os.fstat( fa_genomes.fileno() )
				if filesize.st_size > split_size:
					fa_genomes.close()
					file_count_g += 1
					fa_genomes = open(genomes_string.format(file_count_g), 'wt')

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
