'''
Copy giant Martin multifasta into directory structure needed for CLARK
Run from desired destination directory (/clark_genomes/Custom)
'''

import os
import sys
import lanthpy

def mkdir(line):
	''' Make directory for species in clark folder '''
	dirname = lanthpy.single_name_cleanup(line.partition('|')[0]) # this names the directory by strain
	clark_dir = '.'

	clark = os.path.expanduser(os.path.join(clark_dir,dirname))

	if not os.path.exists(clark):
		os.makedirs(clark)

	return clark

def mkfile(filename,dirname):
	new_filename = os.path.join(dirname, lanthpy.single_name_cleanup(filename.split('|')[-1].strip('_')) +".fna") # this names the file by strain
	#new_filename = os.path.join(dirname, lanthpy.single_name_cleanup(filename.partition('gi|')[2].partition('|')[0] +".fna") # this names the file by GI

	single_fasta = open(new_filename, 'wt')

	# print as a list for tax_id lookup
	print new_filename
	return single_fasta

def main():
	source_file = ' '.join(sys.argv[1:])
	source_file.strip(" ")

	mfa = open(source_file,'r')
	for line in mfa:
		if line[0] == '>': # fasta name
			fdir = mkdir(line)
			ffile = mkfile(line,fdir)
			ffile.write('>gi|' + line.partition('gi|')[2]) # remove the strain label at the front of the header because clark is stupid
		else:
			ffile.write(line)


if __name__ == '__main__':
	main()
