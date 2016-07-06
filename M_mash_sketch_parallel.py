# Build mash index in parallel

import subprocess
import string

for letter in string.ascii_uppercase+"_":
	print letter
	subprocess.call("mash sketch -o ensembl_bact_{0} {0}*.dna.genome.fa -s 100000 -k 21 &".format(letter), shell=True)
