# Run bowtie2 and express on every fastq in a directory, save results with a suffix

import sys
import os
import subprocess

suffix = ""

# get list of fastqs in directory
file_list = [f for f in os.listdir('.') if f.endswith('.fastq')]

# run bowtie2 on each fastq
for f in file_list:
	basename = f[:-6]
	print "Running bowtie2 on {0}".format(basename)
#	proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-k 1000', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x actual_fused_refseq_mrna', '-U ' + f, '-S' + basename + suffix + '.SAM'])
	proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x fused_ensembl_cdna', '-U ' + f, '-S' + basename + suffix + '.SAM'])
	proc.wait()

# run express on each SAM file, in parallel
for f in file_list:
	basename = f[:-6]
	print "Running express on {0}".format(basename)

#	proc = subprocess.Popen(['express', '-o' + basename + suffix, 'actual_fused_refseq_mrna.fa', basename + suffix + '.SAM'])
	proc = subprocess.Popen(['express', '-o' + basename + suffix, 'fused_ensembl_cdna.fa', basename + suffix + '.SAM'])
