# Move all files matching species in text file up a directory

import csv
import os
import pprint

with open('species_hits.txt','r') as hits_file:
	hits_csv = csv.reader(hits_file, delimiter=',')
	hits = [r[0] for r in hits_csv if float(r[1])>1000]

failed_hits = []
for sp in hits:
	cleanup = string.maketrans('-()+/','_____')
	names = [sp.replace(" ","_").lower(), sp.replace(" ","_"), sp.replace(" ","_").lower().translate(cleanup,'.[]'), sp.replace(" ","_").translate(cleanup, '.[]')]
	genuses = [sp.replace(" ","_").lower().partition('_')[0], sp.replace(" ","_").partition('_')[0]]
	success = False
	for name in names: # iterate over possible name forms until one gets a hit
		output = os.system("mv {}* ../".format(name))
		if output == 0:
			success = True
			break
	if not success:
		failed_hits.append(sp)
		for genus in genuses: # iterate over just genus name
			output = os.system("mv {}* ../".format(genus))

pprint.pprint(failed_hits)
print len(failed_hits)
print len(hits)
