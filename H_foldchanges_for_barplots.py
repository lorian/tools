# Create a series of data for bar plots, one for each gene, with the heights pulled from organoid-relative fold change estimates

import sys
import csv
import numpy as np
import os
import string as st

dictlist = {}

# Turn a CSV into a dictionary, indexed by the second column
def make_dict():
	# Turn each of the list of files into a dictionary, stored in a dictionary indexed by the unique part of the filename
	global dictlist
	filelist = ['organoid_control_all_genes.csv','organoid_hES_all_genes.csv','organoid_background_all_genes.csv','organoid_differentiated_all_genes.csv','organoid_teratoma_all_genes.csv','organoid_duodenum_all_genes.csv','organoid_ileum_all_genes.csv','organoid_rectum_all_genes.csv']

	for fn in filelist:
		with open(os.path.join('gene_tables_w21_genenames', fn), 'r') as f:
			data = {rec[1]:rec for rec in csv.reader(f, delimiter=',')}

		dictlist[st.split(fn,'_')[1]] = data

def processgenes(listname,genelist):
	global dictlist
	# Get log2foldchange for each gene, for each sample type
	plotdata = []
	for gene in genelist:
		colnames = []
		genedata = []
		genedata.append(gene)
#		print gene
		for datakey,data in dictlist.iteritems():
#			print "{0}: {1}".format(datakey,data[gene][6])
			colnames.append(datakey)
			try:
				genedata.append(float(data[gene][6]))
			except:
				print data[gene][6]
				genedata.append(data[gene][6])
		plotdata.append(genedata)

	with open('bargraph_data_' +listname+ '.csv', 'wb') as csvfile:
		plotfile = csv.writer(csvfile, dialect='excel')
		plotfile.writerow(['Gene'] + colnames)
		for gd in plotdata:
			plotfile.writerow(gd)

make_dict()
processgenes('aquaporins',["AQP1","AQP2","AQP3","AQP4","AQP5","AQP6","AQP7","AQP8","AQP9","AQP10","AQP11"])
processgenes('maturity',['AXIN2','CDH1','CDH17','CDX2','DEFA5','DEF6','DUOX2','ELF3','EPCAM','KLF4','KLF5','LGR5','LYZ','MUC13','MUC2','SOX9','VIL1','SI','MMP7','PLA2G2A','WNT3A','ATOH1','TFF3','TERT'])
processgenes('adultstem',['LGR5','SOX9','KLF4','KLF5','OLFM4','TERT'])
processgenes('pattern',['BARX1','SFRP1','HHEX','SOX2','DKK1','FABP2','CDX1','CDX2','HOXC5'])
