# Transpose a csv file

import sys
import csv
import numpy as np

filename = sys.argv[1]

with open(filename, 'rb') as csvfile:
	reader = csv.reader( csvfile, 'excel-tab')
	data = []
	for row in reader:
		if len(data) == 0:
			data.append( ['title'] + row )
		else:
			data.append( row )

print data[0][:10]
print data[1][:10]

np_data = np.array( data )
np_data_T = np.transpose(np_data)
print np_data.dtype
np.savetxt('transposed_edit.txt', np_data_T, delimiter='\t', fmt="%s")
