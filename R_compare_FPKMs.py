# Compare FPKMs and est_counts between two results.xprs files

import sys
import csv
import numpy as np

np.set_printoptions(precision=4)

original_filename = sys.argv[1]
new_filename = sys.argv[2]

original_file = open(original_filename,'rU')
original = csv.reader(original_file, 'excel-tab')
original_data = [r for r in original]
original_data = original_data[1:] #remove header row
original_data.sort(key=lambda x:x[1]) # sort by transcript
transcripts = zip(*original_data)[1]

o_est_counts = np.array([float(i) for i in zip(*original_data)[6]])
o_fpkm = np.array([float(i) for i in zip(*original_data)[10]])

new_file = open(new_filename,'rU')
new = csv.reader(new_file, 'excel-tab')
new_data = [r for r in new]
new_data = new_data[1:] #remove header row
new_data.sort(key=lambda x:x[1]) # sort by transcript
transcripts = zip(*new_data)[1]

n_eff_lengths = np.array([float(i) for i in zip(*new_data)[3]])
n_est_counts = np.array([float(i) for i in zip(*new_data)[6]])
n_fpkm = np.array([float(i) for i in zip(*new_data)[10]])

print "ESTIMATED COUNT CHANGES:"
for i in range(len(transcripts)):
	if o_est_counts[i] != n_est_counts[i]:
		print "{0}:\t{1}\t{2}".format(transcripts[i],o_est_counts[i],n_est_counts[i])
print "FPKM CHANGES:"
for i in range(len(transcripts)):
	print "{0}:\tnew FPKM:{1}\tcalculated FPKM:{2}".format(transcripts[i],n_fpkm[i],pow(10,9)*float(n_est_counts[i])/(sum(n_est_counts)*float(n_eff_lengths[i])))

print "TOTAL est_COUNTS:"
print "Original data: {0}".format(sum(o_est_counts))
print "New data: {0}".format(sum(n_est_counts))

