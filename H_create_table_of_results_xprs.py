# Creates a table of a particular column from multiple results.xprs files

import sys
import csv
import numpy as np
import os

count_list = []
fpkm_list = []
try:
	suffix = sys.argv[1] # _quick currently
except:
	suffix = ""

if suffix != "":
	dir_list = [x for x in os.walk('.').next()[1] if x.endswith(suffix)]
else:
	dir_list = [x for x in os.walk('.').next()[1]]


for d in dir_list:
	if "results.xprs" in os.listdir(d) and d != "HC2_duodenum_a":
		if suffix != "":
			basename = d[:-len(suffix)]
		else:
			basename = d
		print "Processing {0}".format(basename)

		results_file = open(os.path.join(d, 'results.xprs'), 'r')
		results = csv.reader(results_file, 'excel-tab')
		results_data = [r for r in results]
		results_data = results_data[1:]
		results_data.sort(key=lambda x:x[1])
		eff_counts = [int(round(float(i if not np.isnan(float(i)) else 0.0))) for i in zip(*results_data)[7]] # for making R input table
		fpkm = [float(i) for i in zip(*results_data)[10]] # for making fpkm master table

		# add sample name to BEGINNING of list
		eff_counts.reverse()
		eff_counts.append(basename)
		eff_counts.reverse()

		fpkm.reverse()
		fpkm.append(basename)
		fpkm.reverse()

		if count_list == []: # first column
			column_labels = [i for i in zip(*results_data)[1]] # gene names
			column_labels.reverse()
			column_labels.append('genes')
			column_labels.reverse()

			count_list.append(column_labels)

		count_list.append(eff_counts)

		if fpkm_list == []: # first column
			fpkm_list.append(column_labels)
		fpkm_list.append(fpkm)

	else:
		print "ERROR: no results.xprs in directory {0}".format(d)

print "Array headers:"
for r in fpkm_list:
	print "Column {0} has {1} rows".format(r[0],len(r))

count_table = np.vstack(count_list)
np.savetxt('r_table_allalign' +suffix+ '.txt', np.transpose(count_table), delimiter='\t', fmt="%s")

fpkm_table = np.vstack(fpkm_list)
np.savetxt('fpkm_table_allalign' +suffix+ '.txt', np.transpose(fpkm_table), delimiter='\t', fmt="%s")
