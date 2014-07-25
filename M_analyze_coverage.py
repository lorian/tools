import pickle
import pprint
import lanthpy

raw_coverage = pickle.load(open('processed_coverage.pickle','r')) # dictionary is: [transcript] = [bp covered, transcript length, counter of coverage depth
diffs = {k:v for k,v in pickle.load(open('testI2_diffs.pickle','r'))}

species = lanthpy.genome_name_cleanup(raw_coverage.keys())
matching_species = [s.replace("_"," ") for s in species]
coverage = {matching_species[i]:raw_coverage[k] for i,k in enumerate(raw_coverage.keys())}

# diffs is only including things above the mean

for transcript in coverage:
	if transcript in diffs.keys():
		print "{}: {} bases covered of {} total length: {:.2%} coverage".format(transcript,coverage[transcript][0],coverage[transcript][1],float(coverage[transcript][0])/float(coverage[transcript][1]))
		print "\tDifference from true value: {}".format(diffs[transcript])
	#for c in depth_coverage:
	#	print "\t{1}X coverage of {0:.2%} of bases".format(float(depth_coverage[c])/float(samfile.lengths[transcript]),c)

# most of "wrong" results are at 90+ coverage?
