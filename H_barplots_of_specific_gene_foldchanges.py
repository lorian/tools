# Create a series of bar plots, one for each gene, with the heights pulled from organoid-relative fold change estimates

import sys
import csv
import numpy as np
import os
import string as st
import matplotlib
#matplotlib.use('macosx') # change to the correct use for your environment
import matplotlib.pyplot as pyp
import random

def do_plots( plot_array, columns = 3, rows = 3, fig_kw = {}, subplot_kw = {}, title = ''):
	""" each item in the plot_array should be a dictionary with the following keys:
		'width' : width of bars
		'color' : color of bars
		'yvals' : array of yvalues
		'xticklabels' : array of label for each yvalue
		'ylabel' : y axis label

		'title' : title for entire plot
		'subplot_kw' : keyword args for creating subplots
		'fig_kw' : keyword args for creating figures
	"""
	num = len(plot_array)
	if num > columns * rows:
		rows = int(np.ceil(num * 1.0 / columns))

	f = pyp.figure( **fig_kw )
	cur_num = 0
	for row_pos in np.arange( rows ):
		for column_pos in np.arange( columns ):
			if cur_num >= num:
				continue
			ax = f.add_subplot( rows, columns, cur_num+1, **subplot_kw )
			plot_single( ax, plot_array[cur_num] )
			cur_num += 1

	#f.title( title )
	f.show()

def plot_single( axes, plot_data_item ):
	axes.bar( np.arange(len(plot_data_item['yvals'])), plot_data_item['yvals'], width = plot_data_item['width'], color=plot_data_item['color'] )
	if 'ylabel' in plot_data_item:
		axes.set_ylabel( plot_data_item['ylabel'] )
	if 'xticklabels' in plot_data_item:
		axes.set_xticklabels( plot_data_item['xticklabels'])
	axes.set_xticks( np.arange( len(plot_data_item['yvals'])) + plot_data_item['width']/2 )

def dummy_plot_array( length ):
	widths = [.2, .5, .8, .95]
	colors = ['r', 'g', 'b']
	yvals_length = [3, 7, 9]
	ylabels = ['cookies', 'length (m)', 'carrots']

	my_array = []
	for i in range(length):
		my_yval_length = random.choice( yvals_length )
		my_yvals = np.random.rand( my_yval_length )
		my_xticklabels = [str(x)[0] for x in my_yvals]
		my_array.append(  { 'yvals': my_yvals,
							'xticklabels': my_xticklabels,
							'ylabel': random.choice( ylabels ),
							'color': random.choice( colors ),
							'width': random.choice( widths ) } )
	return my_array

# recommend calling subplot_adjust to edit subplot parameters:
# defaults:
# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9	# the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9	  # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots

# Turn a CSV into a dictionary, indexed by the second column
def make_dict(filename):
	with open(os.path.join('gene_tables_w21_genenames', filename), 'r') as f:
		data = {rec[1]:rec for rec in csv.reader(f, delimiter=',')}
		return data

def processgenes(genelist):
	# Turn each of the list of files into a dictionary, stored in a dictionary indexed by the unique part of the filename
	filelist = ['organoid_control_all_genes.csv','organoid_hES_all_genes.csv','organoid_background_all_genes.csv','organoid_differentiated_all_genes.csv','organoid_teratoma_all_genes.csv']
	dictlist = {}
	for f in filelist:
		data = make_dict(f)
		dictlist[st.split(f,'_')[1]] = data
	print dictlist.keys()

	# Get log2foldchange for each gene, for each sample type
	plotdata = []

	for gene in genelist:
		datax = []
		datay = []
		for datakey,data in dictlist.iteritems():
			print "{0}: {1}".format(datakey,data[gene][6])
			datax.append(datakey)
			datay.append(float(data[gene][6]))
		plotdata.append({ 'yvals': datay, 'xticklabels': datax, 'title': gene, 'width': 0.8, 'color': 'black'})

		do_plots(plotdata, columns = 3, fig_kw = {'figsize':(6,4)})

def do_example():
	d = dummy_plot_array(12)
	do_plots( d, columns = 5, fig_kw = {'figsize':(10,6) } )  # 10 in by 6 in figure

genelist = ["NES","NEUROG3","PAX6","SOX1","SOX10"]
processgenes(genelist)
