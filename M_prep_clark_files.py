'''
Copy multifastas into directory structure needed for CLARK
'''

import os
import sys
import lanthpy

def mkdir(dirname):
	''' Make directory for species in clark folder '''
	clark = os.path.expanduser(os.path.join('~/scratch/clark/',
				lanthpy.genome_name_cleanup([dirname.partition(".")[0]])[0]))

	if not os.path.exists(clark):
		print "Making directory {0}".format(clark)
		os.makedirs(clark)

	return clark

def mkfile(filename,dirname):
	single_fasta = open(os.path.join(dirname,
			lanthpy.genome_name_cleanup([filename.rstrip().split('|')[-1]])[0] +".fna"), 'wt')
	return single_fasta


def main():
	dirname = ' '.join(sys.argv[1:])
	dirname.lstrip(" ")

	file_list = os.listdir(dirname)
	for f in file_list:
		if f.endswith('fa') or f.endswith('fasta'):
			fdir = mkdir(f)

			mfa = open(os.path.join(dirname,f),'r')
			for line in mfa:
				if line[0] == '>': # fasta name
					ffile = mkfile(line,fdir)
					ffile.write(line)
				else:
					ffile.write(line)


if __name__ == '__main__':
	main()
