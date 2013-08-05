# Change transcript IDs into gene IDs, combine abundances, and recalculate p-values

import sys
import csv
import numpy as np
import os

with open('../ensembl_gene_ids.txt','r') as id_file:
	ids = csv.reader(id_file, 'excel-tab')
	id_data = [r for r in ids]
	id_list = dict(id_data)

#dir_list = [x for x in os.walk('.').next()[1]]
dir_list = ['control_background_all_genes.csv',]
for d in dir_list:
	if d.endswith('.csv'):
		basename = d.rsplit('.')[0]
		print "Processing {0}".format(basename)

		with open(d,'r') as results_file:
			results = csv.reader(results_file, 'excel')
			results_data = [r for r in results]

		updated_results = dict()
		repeated_genes = 0
		new_genes = 0
		header = ""
		for r in results_data:
			try:
				g = id_list[r[1]]
			except:
				g = r[1]
				with open('debug.txt','a') as debug:
					debug.write(str(g)+"\n")

			if g in updated_results:
#				if g not in duplicate_results:
#					duplicate_results[g] = updated_results[g]
				try:
					updated_results[g] = [g] + [float(a)+float(b) for a,b in zip(updated_results[g][1:],r[2:])]
				except:
					pass #lines with 'NA' are all 0 anyway
#				duplicate_results[g] = duplicate_results[g] + r
				repeated_genes = repeated_genes +1
			else:
				if g == 'id':
					header = r[1:]
				else:
					updated_results[g] = [g] + r[2:]
					new_genes = new_genes +1

		print "Unique genes: {0}\nRepeated genes: {1}".format(new_genes,repeated_genes)

		with open(basename+ '_genes.csv','wb') as csvfile:
			update_writer = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
			update_writer.writerow(header)
			for r in updated_results:
				update_writer.writerow(updated_results[r])

