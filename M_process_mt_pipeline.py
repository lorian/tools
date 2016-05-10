# Move all transcriptomes matching species in text file up a directory

import csv
import os
import pprint

with open('species_hits.txt','r') as hits_file:
	hits_csv = csv.reader(hits_file, delimiter=',')
	hits = [r[0] for r in hits_csv if float(r[1])>100]

failed_hits = []
for sp in hits:
	name = sp.lower().replace(" ","_")
	output = os.system("ls {}*".format(name))
	if output != 0: # ls can't match
		name = name.replace('.',"").replace("-","_").replace("(","_").replace(")","_").replace("+","_").replace('[',"").replace(']',"").replace("/","_")
		output = os.system("ls {}*".format(name))
		if output != 0:
			name = sp.partition('_')[0] # get whole genus
			output = os.system("ls {}*".format(name))
			failed_hits.append(sp)

pprint.pprint(failed_hits)
print len(failed_hits)
