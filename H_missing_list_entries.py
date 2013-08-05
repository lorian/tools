# Detect missing transcripts in a list

import sys
import csv
import numpy as np
import os
from operator import itemgetter, attrgetter

with open('ensembl_gene_ids.txt','r') as id_file:
	ids = csv.reader(id_file, 'excel-tab')
	id_data = [r for r in ids]
	id_data.sort(key=itemgetter(0))
	id_list = [i for i in zip(*id_data)[0]] # second column is gene IDs

id_dict = dict(zip(id_list,[0 for x in range(len(id_list))]))

d = 'fpkm_table_ensembl.txt'
basename = d.rsplit('.')[0]
print "Processing {0}".format(basename)

with open(d,'r') as results_file:
	results = csv.reader(results_file, 'excel-tab')
	results_data = [r for r in results]

mfile = open('debug_missing.txt','w')
pfile = open('debug_present.txt','w')
for r in results_data:
	if not r[0] in id_dict:
		mfile.write(str(r[0])+'\n')
	else:
		pfile.write(str(r[0])+'\n')
