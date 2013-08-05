# Mahalanobis distance

import sys
import csv
import numpy as np

comp_file = open('r_PCA_components_no21.txt', 'r')
comps = csv.reader(comp_file, 'excel-tab')
comp_data = np.array([r for r in comps]) # column 0: transcript IDs, col 1: component 1, col 2: component 2

expr_file = open('r_ensembl_vsd_no21.txt', 'r')
exprs = csv.reader(expr_file, 'excel-tab')
expr_data = [r for r in exprs] # title is missing for gene column 1, so all column titles are off by one


# Samples to find distance between
sA = 'DH5_differentiated_line2'
sB = 'DH4_organoid38.23_line2'

expr_A = [float(i) for i in zip(*expr_data[1:])[expr_data[0].index(sA) +1]] # to account for missing column title
expr_B = [float(i) for i in zip(*expr_data[1:])[expr_data[0].index(sB) +1]] # to account for missing column title

# Calculate distance between every gene of the two samples on each component
c1_diff = 0
c2_diff = 0
for i in range(len(expr_A)):
	diff = expr_A[i] - expr_B[i]
	c1_diff_i = abs(diff*float(comp_data[i][1]))
#	print "({0} - {1})*{2} = {3}".format(expr_A[i],expr_B[i],comp_data[i][1],c1_diff_i)
	c1_diff = c1_diff + c1_diff_i

	c2_diff_i = abs(diff*float(comp_data[i][2]))
	c2_diff = c2_diff + c2_diff_i

print c1_diff
print c2_diff

# Weighted average them together to get single distance number
comp1_var = 0.8655761
comp2_var = 0.04134229

distance = comp1_var*c1_diff + comp2_var*c2_diff
print distance
