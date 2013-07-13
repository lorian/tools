# Create a bowtie2 index for every .mfa file in a directory

import sys
import os
import subprocess

dirname = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if dirname == "":
		dirname = arg #avoid extra space at beginning
	else:
		dirname = dirname + " " + arg

file_list = os.listdir(dirname)

for f in file_list:
	if f.find('.mfa') != -1:
		basename = f[:-4]
		print "Making index for {0}".format(f)

		proc = subprocess.Popen(['bowtie2-build', f, basename],stdout=subprocess.PIPE)
		output = proc.communicate()
		if proc.poll() is None:
			raise IOError("bowtie2-build process did not terminate.")
