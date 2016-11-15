import csv
import argparse
import cPickle
import os

parser = argparse.ArgumentParser(description='Look up KEGG pathways from kallisto abundance estimations of genes')
parser.add_argument('filename', help='kallisto abundance.tsv file')
args = parser.parse_args()

if os.path.exists('gene_annotation.pickle'):
	print "Loading annotation dict..."
	annotation_dict = cPickle.load(open('gene_annotation.pickle','rb'))
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
	ortho_files = [f for f in os.listdir('.') if f.startswith("meta.ortho_")]
	ortho_files = ['test.ortho',]
	annotation_dict = dict()
	for f in ortho_files:
		with open(f,'r') as anno_file:
			anno_csv = csv.reader(anno_file, delimiter=' ', skipinitialspace=True)
			anno_data = [r for r in anno_csv]
			anno_dict = {line[1]:pathway_dict[line[5]] for line in anno_data if line[5] in pathway_dict.keys()} # gene id = pathway category
			cPickle.dump(anno_dict,open(f+'.pickle','wb'))
		annotation_dict.update(anno_dict)
		print "Processed {}".format(f)
	print "Saving annotation dict..."
	cPickle.dump(annotation_dict,open('gene_annotation.pickle','wb'))

# save dict as file
with open('gene_annotation.txt','w') as dict_file:
	for k,v in annotation_dict.items():
		dict_file.write(k +'|'+ v +'\n')
		
'''
# get gene IDs from kallisto output
input_file = open(args.filename,'r')
input_csv = csv.reader(input_file, 'excel-tab')
input_data = [r for r in input_csv] #[0] is gene ID, [3] is counts
genes = dict()
genes = {r[0].partition(':')[2].partition('_')[0]:float(r[3]) for r in input_data} # gene:count

for g,c in iteritems(genes):
	try:
		category = annotation_dict[g]

anno_genes = [annotation_dict[g] for g in genes.keys() if g in annotation_dict.keys()]
print anno_genes
'''

# annotate kallisto output
anno_kallisto = open('annotated_'+args.filename,'w')
with open(args.filename,'r') as input_file:
	for line in input_file:
		gene = line.partition(':')[2].partition('_')[0]
		try:
			anno_kallisto.write(annotation_dict[gene] +'|'+ line)
		except:
			anno_kallisto.write(line)
		
anno_kallisto.close()
	
