# Run against mash index in parallel

import subprocess
import string

for letter in string.ascii_uppercase+"_":
	print letter
	subprocess.call('time mash dist ensembl_bact_{0}.msh ~/scratch/illumina_100species_trimmed.all.fq > mash_i100_3_{0}.txt &'.format(letter), shell=True)
