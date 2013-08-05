# Run bowtie2 and express on every fastq in a directory, save results with a suffix

import sys
import os
import subprocess
import tempfile

suffix = ""

# get list of fastqs in directory
file_list = [f for f in os.listdir('.') if f.endswith('.SAM')]

'''
# run bowtie2 on each fastq
for f in file_list:
	basename = f.rsplit('.',1)[0]
	print "Running bowtie2 on {0}".format(basename)
	#proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-k 1000', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x actual_fused_refseq_mrna', '-U ' + f, '-S' + basename + suffix + '.SAM'])
	proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x fused_ensembl_updated_cdna', '-U ' + f, '-S' + basename + suffix + '_newfused.SAM'])
	proc.wait()

# run express on each SAM file, in parallel
for f in file_list:
	basename = f.rsplit('_',1)[0]
	print "Running express on {0}".format(basename)

	#proc = subprocess.Popen(['express', '-o' + basename + suffix, 'actual_fused_refseq_mrna.fa', basename + suffix + '.SAM'])
	proc = subprocess.Popen(['express', '-f0.75', '-B20', '-o' + basename + suffix, 'fused_ensembl_updated_cdna.fa', basename + suffix + '_newfused.bam'])

# run rexpress on each file, in parallel
'''
for f in file_list:
	basename = f.rsplit('.',1)[0]
	if not os.path.isfile(basename+ ".bam") or os.path.getsize(basename+ '.bam') < 10:
		print "Baming file {0}".format(basename)
		os.system('samtools view -bhS ' +basename+ '.SAM > ' +basename+ '.bam')
	if not os.path.isfile(basename+ '_sorted.bam'):
		print "Sorting file {0}".format(basename)
		os.system('samtools sort ' +basename+ '.bam ' +basename+ '_sorted')
	if not os.path.isfile(basename+ '_sorted.bam.bai'):
		print "Indexing file {0}".format(basename)
		os.system('samtools index ' +basename+ '_sorted.bam')

#	print "Running rexpress on {0}".format(basename)
#	proc = subprocess.Popen(['python', 'rexpress.py', basename +'_sorted.bam', basename +'/results.xprs', basename +'/params.xprs', f, 'fused_ensembl_cdna.fa', 'fused_ensembl_updated_cdna.fa', '-g', '-d'+ basename])

