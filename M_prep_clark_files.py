'''
Copy multifastas into directory structure needed for CLARK
'''

import os
import sys
import lanthpy

def mkdir(dirname):
	''' Make directory for species in clark folder '''
	clark_dir = '~/scratch/clark/clark_genomes/Custom/'

	clark = os.path.expanduser(os.path.join(clark_dir,
				lanthpy.single_name_cleanup(dirname.rpartition(".")[0])))

	if not os.path.exists(clark):
		print "Making directory {0}".format(clark)
		os.makedirs(clark)

	return clark

def mkfile(filename,dirname):
	single_fasta = open(os.path.join(dirname,
			lanthpy.single_name_cleanup(filename.rstrip().split('|')[-1]) +".fna"), 'wt')
	return single_fasta

def main():
	source_dir = ' '.join(sys.argv[1:])
	source_dir.lstrip(" ")

	file_list = os.listdir(source_dir)
	for f in file_list:
		if f.endswith('fa') or f.endswith('fasta'):
			print "File {}".format(f)
			fdir = mkdir(f)

			mfa = open(os.path.join(source_dir,f),'r')
			for line in mfa:
				if line[0] == '>': # fasta name
					print "Creating new file for {}".format(line)
					ffile = mkfile(line,fdir)
					ffile.write(line)
				else:
					ffile.write(line)


if __name__ == '__main__':
	main()
