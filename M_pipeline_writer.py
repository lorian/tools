"""
Takes file names and parameters and what's known about existing files and writes a full pipeline shell script
"""
# Todo: check that express directory is empty and delete it
	# make all variables into a parameter file
	# either merge files in a subdirectory or make a variable name to avoid conflicts
	# delete combined_header.txt
	# delete justbody files

import os

use_existing_files = True


def missing_file(filename):
	if use_existing_files and os.path.exists(filename):
		return False
	else:
		return True
		# consider adding lines to the shell script to delete files that will
		# be regenerated


def missing_sam_and_bam(sam_filename):
	bam_filename = insert_suffix(sam_filename, "", 'bam')
	if missing_file(sam_filename) and missing_file(bam_filename):
		return True
	else:
		return False


def insert_suffix(filename,suffix,extension=""):
	""" Insert a string into a filename, and/or change the file extension """
	if extension:
		return "{0}{1}.{2}".format(filename.rpartition('.')[0], suffix,
									extension)
	else:
		return "{0}{1}.{2}".format(filename.rpartition('.')[0], suffix,
									filename.rpartition('.')[2])


def convert_bam_to_sam(sam_filename):
	bam_filename = insert_suffix(sam_filename, "", 'bam')
	if missing_file(sam_filename) and not missing_file(bam_filename):
		return 'samtools view -h {0} > {1} \\\n&& '.format(bam_filename,sam_filename)
	else:
		return False


def main():
	# Filename constants (convert this into a parameter file)
	test_basename = 'testP'
	express_outputname = 'testP1' # for keeping track of different express runs
	version = '2.0'
	cores = 20
	raw_fasta_file = 'illumina_400genomes_named.mfa'
	fastq_file_r1 = "illumina_400species.1.fq.gz"
	fastq_file_r2 = "illumina_400species.2.fq.gz"
	express_cycles = 50
	express_f = 0.90
	nice = False # add in hooks for this

	script = open("M_pipeline_{0}.sh".format(express_outputname),'w')
	script.write("# Version {0} -- used for {1}\n"
				.format(version,express_outputname))

	'''
		# Removes pound sign from fasta file
		cleaned_fasta_file = insert_suffix(raw_fasta_file, '_poundless') # Martin_etal_TextS3_13Dec2011_poundless.fasta
		if missing_file(cleaned_fasta_file):
			script.write('cat {0} | tr -d \\# > {1}\n&& '
						.format(raw_fasta_file,cleaned_fasta_file))

		# Sorts fasta by name of genome
		if missing_file(sorted_fasta_file):
			script.write('sed -e \'/>/s/^/@/\' -e \'/>/s/$/#/\' {0} | tr -d "\\n" | tr "@" "\\n" | sort -t "|" -k1n | tr "#" "\\n" | sed -e \'/^$/d\' > {0} \\\n&& '
						.format(cleaned_fasta_file,sorted_fasta_file))

		# Breaks fasta into bowtie-sized chunks for indexing
		if any(( [missing_file(f) for f in all_fastas] )):
			script.write('python ~/tools/M_process_metagenome_reference.py {0} \\\n&& '
						.format(sorted_fasta_file))
	'''

	# Builds bowtie2 index
	if any((missing_file(insert_suffix(raw_fasta_file, "", "1.bt2")),
			missing_file(insert_suffix(raw_fasta_file, "", "2.bt2")),
			missing_file(insert_suffix(raw_fasta_file, "", "3.bt2")),
			missing_file(insert_suffix(raw_fasta_file, "", "4.bt2")),
			missing_file(insert_suffix(raw_fasta_file, "", "rev.1.bt2")),
			missing_file(insert_suffix(raw_fasta_file, "", "rev.2.bt2")))):
		script.write('bowtie2-build {0} {1} \\\n&& '
					.format(raw_fasta_file, raw_fasta_file.rpartition(".")[0]))

	# Run bowtie2
	sam_file = test_basename + ".SAM"
	if missing_sam_and_bam(sam_file):
		script.write('bowtie2 -a -t -p {0} --local -x {1} -1 {2} -2 {3} -S {4} \\\n&& '
					.format(cores, raw_fasta_file.rpartition(".")[0], fastq_file_r1,
					fastq_file_r2, sam_file))

	# Convert to bam
	bam_file = insert_suffix(sam_file,"","bam") # testL_J.bam
	if missing_file(bam_file):
		script.write('samtools view -bhS {0} > {1} \\\n&& '.format(sam_file,bam_file))

	# Sort
	sorted_bam_file = insert_suffix(bam_file, "_sorted") # testL_J_sorted.bam
	if missing_file(sorted_bam_file):
		script.write('ulimit -n 5000 \\\n&& ')
		script.write('samtools sort -n {0} {1} \\\n&& '
					.format(bam_file, sorted_bam_file.rpartition(".")[0]))

	# Run express (assuming this will always be run)
	script.write('express -f {0} -o {1} --max-indel-size 100 -B {2} {3} {4} \\\n&& ' #apparently the space after -o works now
				.format(express_f, express_outputname, express_cycles, raw_fasta_file, sorted_bam_file))

	# Rename express output and move to parent directory
	script.write('mv {0}/results.xprs {0}_results.xprs \\\n&& '
				.format(express_outputname)) # testL2_results.xprs
	script.write('mv {0}/params.xprs {0}_params.xprs'
				.format(express_outputname)) # testL2_params.xprs


if __name__ == '__main__':
	main()

