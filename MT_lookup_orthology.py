import csv
import argparse
import cPickle
import os

parser = argparse.ArgumentParser(description='Look up KEGG pathways from kallisto abundance estimations of genes')
parser.add_argument('filename', help='kallisto abundance.tsv file')
args = parser.parse_args()

if os.path.exists('gene_annotation.pickle'):
	print "Loading annotation dict..."
	annotation_dict = cPickle.load(open('gene_annotation.pickle'),'rb')
else:
	# get dictionary of KO ids to second-tier pathway names
	pathway_file = open('ko00001.keg','r')
	pathway_csv = csv.reader(pathway_file, delimiter=' ', skipinitialspace=True)
	pathway_data = [r for r in pathway_csv]
	pathway_data = pathway_data[8:] # remove header lines

	pathway_dict = dict()
	for line in pathway_data:
		if line[0] == 'B':
			#B  <b>Cell motility</b>
			pathway = ' '.join(line).partition('<b>')[2].partition('</b>')[0]
		elif line[0] == 'D':
			pathway_dict[line[1]] = pathway


	# get dict of gene ids to KO ids
	annotation_file = open('meta.ortho2','r')
	annotation_csv = csv.reader(annotation_file, delimiter=' ', skipinitialspace=True)
	annotation_data = [r for r in annotation_csv]
	annotation_dict = {line[1]:pathway_dict[line[5]] for line in annotation_data if line[5] in pathway_dict.keys()} # gene id = pathway category

	print "Saving annotation dict..."
	cPickle.dump(annotation_dict,open('gene_annotation.pickle','wb'))


# get gene IDs from kallisto output
input_file = open(args.filename,'r')
input_csv = csv.reader(input_file, 'excel-tab')
input_data = [r for r in input_csv]
genes = [g.partition(':')[2].partition('_')[0] for g in zip(*input_data)[0]]


anno_genes = [annotation_dict[g] for g in genes if g in annotation_dict.keys()]
print anno_genes
