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
if dirname == "" or dirname == " ":
	dirname = os.getcwd() #current working directory

file_list = os.listdir(dirname)

for f in file_list:
	if f.find('.mfa') != -1 and not os.path.isfile(os.path.join(dirname,f.rsplit('.', 1)[0]+ ".1.bt2")):
		basename = f.rsplit('.', 1)[0]

		print "Making index for {0}".format(f)

		proc = subprocess.Popen(['bowtie2-build', os.path.join(dirname,f), os.path.join(dirname,basename)],stdout=subprocess.PIPE)
		output = proc.communicate()
		if proc.poll() is None:
			raise IOError("bowtie2-build process did not terminate.")
