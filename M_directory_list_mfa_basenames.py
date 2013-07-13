# Creates text file containing base names of .mfa files in a directory

import sys
import os

dirname = ""
iterarg = iter(sys.argv)
next(iterarg) #skip name of function
for arg in iterarg:
	if dirname == "":
		dirname = arg #avoid extra space at beginning
	else:
		dirname = dirname + " " + arg

file_list = os.listdir(dirname)

names = open(dirname + "_names.txt","w")
basenames = []
for f in file_list:
	if f.find('.mfa') != -1:
		basename = f[:-4]
		basenames.append(basename)
		names.write(basename + '\n')

print "{0} file names: {1}".format(len(basenames),basenames)
names.close()
