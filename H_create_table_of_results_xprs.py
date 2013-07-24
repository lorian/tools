# Creates a table of a particular column from multiple results.xprs files

import sys
import csv
import numpy as np
import os

r_list = []
suffix = sys.argv[1] # _ensembl normally
dir_list = [x for x in os.walk('.').next()[1] if x.endswith(suffix)]

for d in dir_list:
	if "results.xprs" in os.listdir(d):
		basename = d[:-len(suffix)]
		print "Processing {0}".format(basename)

		results_file = open(os.path.join(d, 'results.xprs'), 'r')
		results = csv.reader(results_file, 'excel-tab')
		results_data = [r for r in results]
		results_data = results_data[1:]
		results_data.sort(key=lambda x:x[1])
		eff_counts = [int(round(float(i))) for i in zip(*results_data)[7]] # for making R input table


		# add sample name to BEGINNING of list
		eff_counts.reverse()
		eff_counts.append(basename)
		eff_counts.reverse()

		if r_list == []: # first column
			column_labels = [i for i in zip(*results_data)[1]] # gene names
			column_labels.reverse()
			column_labels.append('genes')
			column_labels.reverse()

			r_list.append(column_labels)

		r_list.append(eff_counts)

	else:
		print "ERROR: no results.xprs in directory {0}".format(d)

print "Array headers:"
for r in r_list:
	print "Column {0} has {1} rows".format(r[0],len(r))

r_table = np.vstack(r_list)
np.savetxt('r_table' +suffix+ '.txt', np.transpose(r_table), delimiter='\t', fmt="%s")
