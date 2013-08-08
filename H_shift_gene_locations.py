# Shift gene locations in GFF3 file due to insertion

import sys
import csv
import numpy as np
import os
from operator import itemgetter, attrgetter

filename = 'lgr5_allinserts.gff3'
startloc = int('71833550')
insertsize = int('-71833550')

with open(filename) as gff_file:
	gffs = csv.reader(gff_file, 'excel-tab')
	gff_data = [r for r in gffs]

with open(filename.rsplit('.')[0]+ '_shifted.gff3','wb') as csvfile:
	shifted = csv.writer(csvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
	for r in gff_data:
		new_r = r
		if r and not r[0].startswith('#'):
			if int(r[3])>=startloc:
				new_r[3] = int(r[3]) + insertsize
			if int(r[4])>startloc:
				new_r[4] = int(r[4]) + insertsize
		shifted.writerow(new_r)

#12	Ensembl	exon	71834586	71837748	.	+	2	Name=NtermGFP;Parent=ENST00000266674_Nterm_Ntype

