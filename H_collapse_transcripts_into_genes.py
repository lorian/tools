# Change transcript IDs into gene IDs and combine abundances

import sys
import csv
import numpy as np
import os
from operator import itemgetter, attrgetter

with open('ensembl_gene_ids.txt','r') as id_file:
	ids = csv.reader(id_file, 'excel-tab')
	id_data = [r for r in ids]
	id_list = dict(id_data)

with open('ensembl_gene_names.txt','r') as name_file:
	names = csv.reader(name_file, 'excel-tab')
	name_data = [r for r in names]
	name_list = dict(name_data)

# Replace ensembl ids with gene names, where they exist
named_dict = dict()
for transcript in id_list:
	try:
		name = name_list[id_list[transcript]]
	except:
		named_dict[transcript] = transcript
	if name == "":
		named_dict[transcript] = id_list[transcript]
	else:
		named_dict[transcript] = name

#d = 'fpkm_table_newfused.txt'
d = 'r_table_allalign_a.txt'
basename = d.rsplit('.')[0]
print "Processing {0}".format(basename)

with open(d,'r') as results_file:
	results = csv.reader(results_file, 'excel-tab')
	results_data = [r for r in results]

updated_results = dict()
repeated_genes = 0
new_genes = 0
header = ""
debug = open('debug.txt','w')
for r in results_data:
	try:
		g = named_dict[r[0]]
	except:
		g = r[0]
		debug.write(str(g)+"\n")

	if g in updated_results:
		updated_results[g] = [g] + [float(a)+float(b) for a,b in zip(updated_results[g][1:],r[1:])]
		repeated_genes = repeated_genes +1
	else:
		if g == 'genes':
			header = r
		else:
			updated_results[g] = [g] + r[1:]
			new_genes = new_genes +1

print "Unique genes: {0}\nRepeated genes: {1}".format(new_genes,repeated_genes)

with open(basename+ '_genenames.txt','wb') as csvfile:
	update_writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
	update_writer.writerow(header)
	for r in updated_results:
		update_writer.writerow(updated_results[r])

