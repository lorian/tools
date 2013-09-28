# Merge multiple SAM files together

import sys
import os
import shutil
import pysam

filename = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
filenames = []
for arg in iterarg:
	filenames.append(arg)

samfiles = []
for name in filenames:
	print name
	samfiles.append(pysam.Samfile(name,'r'))

full_references = []
full_lengths = []
for f in samfiles:
	full_references = full_references + [r for r in f.references]
	full_lengths = full_lengths + [r for r in f.lengths]

header_sam = pysam.Samfile('combined_header.sam', mode='wh', referencenames=full_references, referencelengths=full_lengths)
header_sam.close()

full_sam = open('combined_file.sam','w')
shutil.copyfileobj(open('combined_header.sam','r'), full_sam)
for name in filenames:
	os.system('samtools view -S -o {0}_justbody.txt {0}'.format(name))
	shutil.copyfileobj(open(name + '_justbody.txt','r'), full_sam)
full_sam.close()
