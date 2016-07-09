# Build mash index in parallel

import subprocess
import string

for letter in string.ascii_uppercase+"_":
	print letter
	subprocess.call("time mash sketch -o ensembl_bact_{0} {0}*.dna.genome.fa -s 10000 -k 25 2> log_bact_{0} &".format(letter), shell=True)
