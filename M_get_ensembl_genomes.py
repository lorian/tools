# Download all cdna files from ensembl bacteria listings

import os

def get_listing(ftp):
	os.system("wget -O index.html {}".format(ftp))
	with open('index.html','r') as f:
		# get listing of directories from html format
		dir_names = [line.partition("<a href=\"")[2].partition("\">")[0] for line in f]
		dir_names = [d for d in dir_names if d] # drop empty lines that didn't have links
		return dir_names

collection_names = get_listing("ftp://ftp.ensemblgenomes.org/pub/current/bacteria/fasta/")

for c in collection_names:
	print "\t\tGetting species for collection {}".format(c.split('/')[-2])
	species_names = get_listing(c)
	for sp in species_names:
		# change "cdna" and ".cdna.all.fa.gz" here if you want to download DNA
		name = sp.split('/')[-2]
		print "\t\t{}".format(name)
		if not os.path.exists("{}.cdna.all.fa.gz".format(name)):
			os.system("wget {}cdna/*.cdna.all.fa.gz -O {}.cdna.all.fa.gz".format(sp,name))
