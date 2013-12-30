# Output R code to produce an analysis set

import itertools

conditions = ['hES','teratoma','organoid','differentiated','control','background','mixed']

for pair in itertools.combinations(conditions,2):
#	print "ensembl_cds_binom = nbinomTest( ensembl_filt, \"{0}\",\"{1}\" )".format(pair[0],pair[1])
#	print "write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file=\"{0}_{1}_filtered_genes.csv\" )".format(pair[0],pair[1])

	print "ensembl_cds_binom = nbinomTest( ensembl_cds_disp, \"{0}\",\"{1}\" )".format(pair[0],pair[1])
	print "write.csv( ensembl_cds_binom[ order(ensembl_cds_binom$pval), ], file=\"decemberdata_{0}_{1}_all_genes.csv\" )".format(pair[0],pair[1])
