# Replace fasta headers with long name versions

import os

sp_dict = dict()
with open('IMG_ID2SPECIES_LIST.txt','r') as sp_map:
	for line in sp_map:
		segments = line.partition(" ")
		sp_dict[segments[0]] = segments[2].replace('[','').replace(' ','_')

fasta = open('ld_7point5mill_sample.fasta','r')
with open('ld_7point5mill_sample_names.fasta','w') as newfile:
	for line in fasta:
		if line.startswith('>'):
			uid = line[1:].partition('_')[0]
			newfile.write('>'+ line.partition('_')[2].strip() +'_'+ sp_dict[uid])
		else:
			newfile.write(line)
