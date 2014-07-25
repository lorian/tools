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
print filenames

samfiles = []
for name in filenames:
	samfiles.append(pysam.Samfile(name,'r'))

full_references = []
full_lengths = []
for f in samfiles:
	full_references = full_references + [r for r in f.references]
	full_lengths = full_lengths + [r for r in f.lengths]

print "Creating combined header..."
header_sam = pysam.Samfile('combined_header.sam', mode='wh', referencenames=full_references, referencelengths=full_lengths)
header_sam.close()

full_sam = open('combined_file.sam','w')
shutil.copyfileobj(open('combined_header.sam','r'), full_sam)
for name in filenames:
	print "Copying {0}".format(name)
	#if not os.path.isfile(name + '_justbody.txt'): #don't recreate body file if it already exists from previous failed runs
	os.system('samtools view -S -o {0}_justbody.txt {0}'.format(name))
	shutil.copyfileobj(open(name + '_justbody.txt','r'), full_sam)
full_sam.close()
