# Create a series of data for bar plots, one for each gene, with the heights pulled from organoid-relative fold change estimates

import sys
import csv
import numpy as np
import os
import string as st
import math

dictlist = {}
varlist = {}
meanlist = {}

# Turn a CSV into a dictionary, indexed by the second column
def make_dict():
	# Turn each of the list of files into a dictionary, stored in a dictionary indexed by the unique part of the filename
	global dictlist
	global varlist
	global meanlist
	filelist = ['organoid_control_new_genes.csv','organoid_hES_new_genes.csv','organoid_background_new_genes.csv','organoid_differentiated_new_genes.csv','organoid_teratoma_new_genes.csv']

	for fn in filelist:
		with open(os.path.join('decemberdata', fn), 'r') as f:
			data = {rec[1]:rec for rec in csv.reader(f, delimiter=',')}

		dictlist[st.split(fn,'_')[1]] = data

	var_file = open('decemberdata/variances.txt','r')
	varlist = {rec[0]:rec[1:] for rec in csv.reader(var_file, delimiter='\t')}
	mean_file = open('decemberdata/means.txt','r')
	meanlist = {rec[0]:rec[1:] for rec in csv.reader(mean_file, delimiter='\t')}

def processgenes(listname,genelist):
	global dictlist
	global varlist
	global meanlist
	# Get log2foldchange for each gene, for each sample type
	plotdata = []
	for gene in genelist:
		colnames = []
		genedata = []
		genedata.append(gene)
#		print gene
		for datakey,data in dictlist.iteritems():
#			print "{0}: {1}".format(datakey,data[gene][6])
			colnames.append(datakey) #order is actually a constant
			try:
				genedata.append(float(data[gene][6]))
			except:
				genedata.append(data[gene][6]) #inf or -inf

		# Calculate variances
		# Assuming order of var and means is: organoid, differentiated, control, hES, background, teratoma
		varpart = [0 if (meanlist[gene][i] == '0') else float(varlist[gene][i])/(float(meanlist[gene][i]))**2 for i,d in enumerate(varlist[gene])]
		vardata = [math.sqrt((1/math.log(2)**2)*(varpart[i] + varpart[0])) for i,d in enumerate(varpart)]
		colnames = colnames + ['differentiated stdev','control stdev','hES stdev','background stdev','teratoma stdev']
		plotdata.append(genedata + vardata[1:])


	with open('bargraph_data_stdev_' +listname+ '.csv', 'wb') as csvfile:
		plotfile = csv.writer(csvfile, dialect='excel')
		plotfile.writerow(['Gene'] + colnames)
		for gd in plotdata:
			plotfile.writerow(gd)

make_dict()
#processgenes('aquaporins',["AQP1","AQP2","AQP3","AQP4","AQP5","AQP6","AQP7","AQP8","AQP9","AQP10","AQP11"])
#processgenes('maturity',['AXIN2','CDH1','CDH17','CDX2','DEFA5','DEF6','DUOX2','ELF3','EPCAM','KLF4','KLF5','LGR5','LYZ','MUC13','MUC2','SOX9','VIL1','SI','MMP7','PLA2G2A','WNT3A','ATOH1','TFF3','TERT'])
#processgenes('adultstem',['LGR5','SOX9','KLF4','KLF5','OLFM4','TERT'])
#processgenes('pattern',['BARX1','SFRP1','HHEX','SOX2','DKK1','FABP2','CDX1','CDX2','HOXC5'])
#processgenes('heatmap',["AFP","ANPEP","ASCL2","AXIN2","BARX1","CA2","CA2","CD34","CDH1","CDH17","CDX1","CDX2","CDX4","CLCA4","CLDN6","CTNNB1","CXCR4","DEF6","DKK1","DUOX2","DUOX2","ELF3","EOMES","EPCAM","FABP1","FABP2","FABP2","FOXA1","FOXA2","GATA4","GATA4","GATA6","GDF1","GDF3","GSC","GUCA2A","HHEX","HNF1B","HNF4A","HOXC5","KIT","KLF4","KLF5","KRT19","KRT20","LGR5","LYZ","MEP1A","MGAM","MIXL1","MKI67","MMP7","MUC13","MUC17","MUC2","NKX2-1","OLFM4","PDX1","PLA2G2A","REG4","SALL4","SCNN1A","SFRP1","SI","SLC15A1","SLC5A1","SOX17","SOX2","SOX7","SOX9","TERT","TFF2","TFF3","VIL1","VILL","WNT3A"])
#processgenes('qPCRcompare',["LGR5","SOX9","CDX2","CDX1","KLF4","KLF5","CDH1","MUC2","VIL1","SI","ALPI","CDH2","NANOG","GATA4","FGF8","SOX17","ACTB","DEFA5","DEFA6"])
processgenes('finalset',["AXIN2","LGR5","ATOH1","MSI1","EPHB2","DCLK1","LRIG1","TERT","AXIN2","OLFM4","ASCL2","MIXL1","NES","TFF3"])
