import pickle
import pprint
import lanthpy

raw_coverage = pickle.load(open('processed_coverage.pickle','r')) # dictionary is: [transcript] = [bp covered, transcript length, counter of coverage depth
diffs = {k:v for k,v in pickle.load(open('testE_L_diffs.pickle','r'))}

species = lanthpy.genome_name_cleanup(raw_coverage.keys())
coverage = {species[i]:raw_coverage[k] for i,k in enumerate(raw_coverage.keys())}
print "diffs transcripts: {}\ncoverage transcripts: {}".format(len(diffs.keys()),len(coverage.keys()))

diff_coverage = {} # dictionary that contains a list of all the coverages of transcripts with that diff
for transcript in diffs.keys():
	if transcript in coverage:
		cov = float(coverage[transcript][0])/float(coverage[transcript][1])
		#print "{}: {} bases covered of {} total length: {:.2%} coverage".format(transcript,coverage[transcript][0],coverage[transcript][1],cov)
		#print "\tDifference from true value: {}".format(diffs[transcript])
		diff_coverage.setdefault(round(diffs[transcript]),[]).append(round(100*cov))
	else:
		print "{} present in diffs, not present in coverage analysis.".format(transcript)
	#for c in depth_coverage:
	#	print "\t{1}X coverage of {0:.2%} of bases".format(float(depth_coverage[c])/float(samfile.lengths[transcript]),c)
'''
for diff in sorted(diff_coverage.keys()):
	num = len(diff_coverage[diff])
	avg = sum(diff_coverage[diff])/num
	print "There were {} diffs of {}, with an average coverage of {}".format(num,diff,avg)
'''

nonzero_coverages = [v for k,v in diff_coverage.items() if k != 0]
print nonzero_coverages
# most of "wrong" results are at 90+ coverage?
# count results in each diff range, see what coverage is like

[[100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0, 100.0, 100.0], [100.0, 100.0], [98.0, 100.0], [100.0], [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 98.0, 98.0, 100.0, 100.0, 100.0], [100.0, 100.0, 100.0, 100.0, 100.0], [100.0, 100.0, 98.0], [100.0], [100.0], [100.0], [100.0, 100.0], [100.0], [100.0], [99.0], [100.0], [100.0], [2.0, 100.0], [100.0], [100.0, 100.0, 100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0], [100.0, 100.0], [97.0], [100.0, 100.0], [1.0], [100.0], [100.0], [73.0, 27.0, 70.0, 97.0, 95.0, 77.0, 70.0, 4.0, 87.0, 94.0, 59.0, 99.0, 89.0, 85.0, 100.0, 92.0, 87.0, 33.0, 4.0, 90.0, 100.0, 58.0, 56.0, 82.0, 74.0, 98.0, 95.0, 96.0, 2.0, 92.0, 71.0, 89.0, 36.0, 88.0, 79.0, 61.0, 55.0, 91.0, 60.0, 100.0, 95.0, 85.0, 85.0, 90.0, 77.0, 45.0, 95.0, 94.0, 0.0, 91.0, 99.0, 97.0, 9.0, 86.0, 65.0, 61.0, 9.0, 84.0, 2.0, 98.0, 3.0, 79.0, 70.0, 52.0, 52.0, 60.0, 60.0, 6.0, 59.0, 80.0, 83.0, 63.0, 97.0, 83.0, 60.0, 92.0, 100.0, 58.0, 2.0, 59.0, 57.0, 2.0, 47.0, 72.0, 95.0, 92.0, 97.0, 7.0, 57.0, 59.0, 88.0, 97.0, 94.0, 61.0, 100.0, 58.0, 60.0, 56.0, 1.0, 78.0, 46.0, 66.0, 66.0, 15.0, 39.0, 92.0, 58.0, 48.0, 68.0, 2.0, 62.0, 74.0, 63.0, 33.0, 98.0, 86.0, 60.0, 92.0, 59.0, 100.0, 63.0, 87.0, 94.0, 60.0, 99.0, 75.0, 10.0, 56.0, 99.0, 100.0, 27.0, 61.0, 54.0, 94.0, 75.0, 12.0, 92.0, 69.0, 97.0, 94.0, 60.0, 5.0, 79.0, 91.0, 61.0, 100.0, 82.0, 80.0, 65.0, 96.0, 66.0, 48.0, 93.0, 57.0, 53.0, 74.0, 57.0, 52.0, 83.0, 13.0, 97.0, 57.0, 52.0, 97.0, 7.0, 59.0, 4.0, 60.0, 60.0, 28.0, 41.0, 2.0, 58.0, 60.0], [100.0, 100.0], [100.0], [100.0, 100.0, 100.0], [100.0], [100.0], [99.0], [99.0, 97.0, 100.0, 96.0], [100.0, 100.0, 100.0], [100.0, 100.0, 100.0], [100.0], [100.0], [100.0], [100.0, 100.0, 100.0], [100.0], [100.0], [100.0], [100.0]]
