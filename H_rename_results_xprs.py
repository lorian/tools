# Rename results.xprs files based on directory, and moves _noa ones to GEO directory

import os
import shutil

dir_list = [x for x in os.walk('.').next()[1]]
for d in dir_list:
	if d == "GEO":
		continue
	if "results.xprs" in os.listdir(d):
		print "Renaming results.xprs in {0} to {0}.xprs".format(d)
		os.renames(os.path.join(d,"results.xprs"),os.path.join(d,d+".xprs"))

	file_list = [x for x in os.walk(d).next()[2]]
	for f in file_list:
		if f.endswith("_noa.xprs"):
			print "Copying {0} to GEO/".format(f)
			shutil.copy(os.path.join(d,f),"GEO")

