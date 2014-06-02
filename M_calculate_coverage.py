# Determine how well-covered by reads and depth each metagenome is
import pysam
import collections

samfile = pysam.Samfile("../metagenome/testK2_short_sorted.bam", "rb" )

transcript_coverage = {}
for pileupcolumn in samfile.pileup():
	transcript_coverage.setdefault(pileupcolumn.tid,[]).append(pileupcolumn.n)
		#print "Coverage for transcript {} of length {} at base {} = {}".format(samfile.getrname(pileupcolumn.tid),samfile.lengths[pileupcolumn.tid],pileupcolumn.pos , pileupcolumn.n)

for transcript in transcript_coverage:
	depth_coverage = collections.Counter(transcript_coverage[transcript])
	print "{}: {} bases covered of {} total length: {:.2%} coverage".format(samfile.getrname(transcript),len(transcript_coverage[transcript]),samfile.lengths[transcript],float(len(transcript_coverage[transcript]))/float(samfile.lengths[transcript]))
	for c in depth_coverage:
		print "\t{1}X coverage of {0:.2%} of bases".format(float(depth_coverage[c])/float(samfile.lengths[transcript]),c)


samfile.close()

