# Transpose a csv file

import sys
import csv
import numpy as np

filename = sys.argv[1]

with open(filename, 'rb') as csvfile:
	reader = csv.reader( csvfile, 'excel-tab')
	data = []
	i = 0
	for row in reader:
		if i == 0:
			first_row = row
		elif i== 1:
			if len(first_row) < len(row): # checks to see if a column title is needed
				data.append(['title'] + first_row)
				data.append(row)
			else:
				data.append(first_row)
				data.append(row)
		else:
			data.append( row )
		i= i+1

np_data = np.array( data )
np_data_T = np.transpose(np_data)
print np_data.dtype
np.savetxt('transposed_edit.txt', np_data_T, delimiter='\t', fmt="%s")
