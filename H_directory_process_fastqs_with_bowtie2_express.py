# Run bowtie2 and express on every fastq in a directory, save results with a suffix

import sys
import os
import subprocess
import tempfile

suffix = "_quick"

# get list of fastqs in directory
file_list = [f for f in os.listdir('.') if f.endswith('.fastq') and not 'contaminated' in f]
#file_list = [f for f in os.listdir('.') if f.endswith('.SAM')]

def bowtie2():
	# run bowtie2 on each fastq
	for f in file_list:
		basename = f.rsplit('.',1)[0]
		if not os.path.isfile(basename + suffix +'.SAM'):
			print "Running bowtie2 on {0}".format(basename)
			#proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-k 1000', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x actual_fused_refseq_mrna', '-U ' + f, '-S' + basename + suffix + '.SAM'])
			proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-a', '-t', '-p 35', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x fused_ensembl_updated_cdna', '-U ' + f, '-S' + basename + suffix +'.SAM'])
#			proc = subprocess.Popen(['../bowtie2-2.1.0/bowtie2', '-a', '-t', '-p 40', '--rdg 6,5', '--rfg 6,5', '--score-min L,-.6,-.4', '-x ENST74_fused', '-U ' + f, '-S' + basename + '_LGR5.SAM'])
			proc.wait()

def express():
	# run express on each SAM file, in parallel
	for f in file_list:
		basename = f.rsplit('.',1)[0]
		print "Running express on {0}".format(basename)
		proc = subprocess.Popen(['express', '-f0.75', '-o' + basename + suffix, 'fused_ensembl_updated_cdna.fa', basename +'_a.SAM'])

def prep_sam():
	# prep each file for alignment viewing, in parallel
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

def filestuff():
	for f in file_list:
		basename = f.split('_',1)[0]
		print "Processing {0}".format(basename)
		os.system('grep -E \'ENST00000266674\' ' +f+ ' > ' +basename+ '.txt')


#bowtie2()
express()
#prep_sam()
