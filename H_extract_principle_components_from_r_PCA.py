# Extract principle components from an r-generated PCA

import sys
import csv
import numpy as np

PCA_file = open('r_PCA_scores_no21.txt', 'r')
PCA = csv.reader(PCA_file, 'excel-tab')

PCA_data = [r for r in PCA]
PCA_data = PCA_data[1:] # Remove first (header) row

comp_data = zip(*PCA_data) # Transpose table
comp_table = np.vstack((comp_data[0],comp_data[1],comp_data[2]))
np.savetxt('r_PCA_components_no21.txt', np.transpose(comp_table), delimiter='\t', fmt="%s") # Save columns 1,2, and 3

'''
# Sorted separate components
PCA_data.sort(key=lambda x:abs(float(x[1])), reverse=True) # Sort based on absolute value of column 2
comp1_data = zip(*PCA_data) # Transpose sorted table
comp1_table = np.vstack((comp1_data[0],comp1_data[1]))
np.savetxt('r_PCA_component1.txt', np.transpose(comp1_table), delimiter='\t', fmt="%s") # Save columns 1 and 2

PCA_data.sort(key=lambda x:abs(float(x[2])), reverse=True) # Sort based on absolute value of column 3
comp2_data = zip(*PCA_data) # Transpose sorted table
comp2_table = np.vstack((comp2_data[0],comp2_data[2]))
np.savetxt('r_PCA_component2.txt', np.transpose(comp2_table), delimiter='\t', fmt="%s") # Save columns 1 and 3
'''
