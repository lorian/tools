# Rename .fna files downloaded from NCBI

import sys
import os
import lanthpy

dirname = lanthpy.get_arg_w_spaces(sys.argv[1:])
subdirs_list = os.listdir(dirname)

for d in subdirs_list:
	file_list = os.listdir(os.path.join(dirname,d))
	print d
	for f in file_list:
		if f.endswith('.fna'):
			with open(os.path.join(dirname,d,f),'r+') as mfa:
				for line in mfa:
					name = '>'+ d.lstrip('_').partition('uid')[0] + line[1:].replace(" ","_")
					mfa.seek(0)
					mfa.write(name)
					break
