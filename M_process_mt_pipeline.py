# Move all transcriptomes matching species in text file up a directory

import csv
import os

with open('species_hits.txt','r') as hits_file:
	hits_csv = csv.reader(hits_file, delimiter=',')
	hits = [r[0] for r in hits_csv if float(r[1])>50]
	print hits
	print len(hits)

for sp in hits:
	os.system("ls {}*".format(sp.lower().replace(" ","_")))
