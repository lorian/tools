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
	test_basename = 'testN'
	express_outputname = 'testN1' # for keeping track of different express runs
	version = '1.4'
	cores = 40
	raw_fasta_file = 'Martin_etal_TextS3_13Dec2011_poundless_fragmented.fasta'
	fastq_file_r1 = "illumina_100species.1.fq.gz"
	fastq_file_r2 = "illumina_100species.2.fq.gz"
	i100_alignment = 'testNi100.SAM'
	i100_fasta = 'illumina_100genomes_fragmented_genomes_0.mfa'
	express_cycles = 20
	express_f = 0.85
	nice = True # add in hooks for this

	script = open("M_pipeline_{0}.sh".format(express_outputname),'w')
	script.write("# Version {0} -- used for {1}\n"
				.format(version,express_outputname))

	# Looks for bowtie-sized chunks for indexing
#	NOTE TEMP CHANGE HERE TO WORK WITH TEST N!
	sorted_fasta_file = insert_suffix(raw_fasta_file, '') # Martin_etal_TextS3_13Dec2011_sorted.fasta
	all_fastas = [insert_suffix(sorted_fasta_file,'_genomes_0','mfa'),
					insert_suffix(sorted_fasta_file,'_genomes_1','mfa'),
					insert_suffix(sorted_fasta_file,'_genomes_2','mfa'),
					insert_suffix(sorted_fasta_file,'_plasmids_0','mfa')] # Martin_etal_TextS3_13Dec2011_sorted_genomes_0.mfa
	if any(( [missing_file(f) for f in all_fastas] )):
		# Go through entire fasta creation process

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

	# Builds bowtie2 indexes
	for fasta in all_fastas: # Martin_etal_TextS3_13Dec2011_sorted_genomes_1.1.bt2
		if any((missing_file(insert_suffix(fasta, "", "1.bt2")),
				missing_file(insert_suffix(fasta, "", "2.bt2")),
				missing_file(insert_suffix(fasta, "", "3.bt2")),
				missing_file(insert_suffix(fasta, "", "4.bt2")),
				missing_file(insert_suffix(fasta, "", "rev.1.bt2")),
				missing_file(insert_suffix(fasta, "", "rev.2.bt2")))):
			script.write('bowtie2-build {0} {1} \\\n&& '
						.format(fasta, fasta.rpartition(".")[0]))

	# Run bowtie2
	all_suffixes = ['_g0', '_g1', '_g2', '_p'] # testL_g0.SAM
	for i,fasta in enumerate(all_fastas):
		if missing_sam_and_bam(test_basename + all_suffixes[i] + ".SAM"):
			script.write('nice -n 19 bowtie2 -a -t -p {0} --local -x {1} -1 {2} -2 {3} -S {4}.SAM \\\n&& '
						.format(cores, fasta.rpartition(".")[0], fastq_file_r1,
						fastq_file_r2, test_basename + all_suffixes[i]))

	# Merge sam files
	merged_sam_file = "{0}_{1}".format(test_basename, i100_alignment.rpartition("test")[2]) # testL_J.SAM
	if missing_sam_and_bam(merged_sam_file):

		# convert any bams into sams for merge
		if convert_bam_to_sam(i100_alignment):
			script.write(convert_bam_to_sam(i100_alignment))
		for i,fasta in enumerate(all_fastas):
			if convert_bam_to_sam(test_basename + all_suffixes[i] + ".SAM"):
				script.write(convert_bam_to_sam(test_basename + all_suffixes[i] + ".SAM"))

		# create i100 alignment if it doesn't exist
		if missing_file(i100_alignment):
			if any((missing_file(insert_suffix(i100_fasta, "", "1.bt2")),
					missing_file(insert_suffix(i100_fasta, "", "2.bt2")),
					missing_file(insert_suffix(i100_fasta, "", "3.bt2")),
					missing_file(insert_suffix(i100_fasta, "", "4.bt2")),
					missing_file(insert_suffix(i100_fasta, "", "rev.1.bt2")),
					missing_file(insert_suffix(i100_fasta, "", "rev.2.bt2")))):
				script.write('bowtie2-build {0} {1} \\\n&& '
						.format(i100_fasta, i100_fasta.rpartition(".")[0]))
			script.write('nice -n 19 bowtie2 -a -t -p {0} --local -x {1} -1 {2} -2 {3} -S {4} \\\n&& '
					.format(cores, i100_fasta.rpartition(".")[0], fastq_file_r1,
					fastq_file_r2, i100_alignment))

		# Merge sams
		script.write('python ~/tools/M_merge_sam_files.py {0} {1} \\\n&& '
					.format(" ".join([test_basename + s + '.SAM' for s in all_suffixes]),i100_alignment))
		script.write('mv combined_file.sam {0} \\\n&& '.format(merged_sam_file))

	# Convert to bam
	bam_file = insert_suffix(merged_sam_file,"","bam") # testL_J.bam
	if missing_file(bam_file):
		script.write('samtools view -bhS {0} > {1} \\\n&& '.format(merged_sam_file,bam_file))

	# Sort
	sorted_bam_file = insert_suffix(bam_file, "_sorted") # testL_J_sorted.bam
	if missing_file(sorted_bam_file):
		script.write('ulimit -n 5000 \\\n&& ')
		script.write('samtools sort -n {0} {1} \\\n&& '
					.format(bam_file, sorted_bam_file.rpartition(".")[0]))

	# Combine all mfa files used
	merged_all_fastas = insert_suffix(sorted_fasta_file, "_i100", "mfa") # Martin_etal_TextS3_13Dec2011_sorted_i100.mfa
	if missing_file(merged_all_fastas):
		script.write('cat {0} {1} > {2} \\\n&& '
					.format(i100_fasta, " ".join([f for f in all_fastas]), merged_all_fastas))

	# Run express (assuming this will always be run)
	script.write('express -f {0} -o {1} --max-indel-size 100 -B {2} {3} {4} \\\n&& ' #apparently the space after -o works now
				.format(express_f, express_outputname, express_cycles, merged_all_fastas, sorted_bam_file))

	# Rename express output and move to parent directory
	script.write('mv {0}/results.xprs {0}_results.xprs \\\n&& '
				.format(express_outputname)) # testL2_results.xprs
	script.write('mv {0}/params.xprs {0}_params.xprs'
				.format(express_outputname)) # testL2_params.xprs


if __name__ == '__main__':
	main()
