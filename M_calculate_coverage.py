# Determine how well-covered by reads and depth each metagenome is
import pysam
import collections
import pickle
import string

samfile = pysam.Samfile("testE_I_sorted.bam", "rb" )
has_pickle = True

if not has_pickle:
	transcript_coverage = {}
	for pileupcolumn in samfile.pileup():
		transcript_coverage.setdefault(pileupcolumn.tid,[]).append(pileupcolumn.n)
			#print "Coverage for transcript {} of length {} at base {} = {}".format(samfile.getrname(pileupcolumn.tid),samfile.lengths[pileupcolumn.tid],pileupcolumn.pos , pileupcolumn.n)

	pickle.dump(transcript_coverage,open('testE_I_pileup.pickle','w'))

else:
	transcript_coverage = pickle.load(open('testE_I_pileup.pickle','r'))

processed_coverage = dict()
for transcript in transcript_coverage:
	depth_coverage = collections.Counter(transcript_coverage[transcript])
	processed_coverage[samfile.getrname(transcript)] = [len(transcript_coverage[transcript]),samfile.lengths[transcript],depth_coverage] # [transcript name] = [covered bases, length, counter]
	#print "{}: {} bases covered of {} total length: {:.2%} coverage".format(samfile.getrname(transcript),len(transcript_coverage[transcript]),samfile.lengths[transcript],float(len(transcript_coverage[transcript]))/float(samfile.lengths[transcript]))
	#for c in depth_coverage:
	#	print "\t{1}X coverage of {0:.2%} of bases".format(float(depth_coverage[c])/float(samfile.lengths[transcript]),c)

pickle.dump(processed_coverage,open('processed_coverage.pickle','w'))

samfile.close()

