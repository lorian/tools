import ftplib

ftp = ftplib.FTP("ftp.ensemblgenomes.org")
ftp.login()
ftp.cwd("pub/current/bacteria/fasta/")
dir_list = [d for d in ftp.nlst() if d.startswith('bacteria')]

for d in dir_list:
	bact_list = ftp.nlst(d)
	for b in bact_list:
		fasta_list = [f for f in ftp.nlst(b + '/dna') if f.endswith('.dna.genome.fa.gz')]
		print fasta_list
		ftp.retrbinary('RETR '+ fasta_list[0], open(fasta_list[0].lstrip(b +'/dna/'), 'wb').write)
