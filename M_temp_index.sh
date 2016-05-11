#!/bin/bash
echo "truth"
cd ~/scratch/ensembl_cdna/truth/
time kallisto index -i new_truth_index *.cdna.all.fa > log_index_new_truth.txt 2>&1
time kallisto quant -i new_truth_index -o kallisto_new_truth --single ~/scratch/mt_sim/ld_7point5mill_sample_processed.fasta -l 200 -s 5 -t 30 > log_align_new_truth.txt 2>&1
mv kallisto_new_truth/abundance.tsv ~/scratch/ensemble_cdna/kallisto_new_truth_abundance.txt
echo "sim test genomes"
indexes=(1 2 3 4 6 7)
FASTQ="ld_7point5mill_sample_processed.fasta"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	time kallisto quant -i singlesp_genome_${INDEX}_index -o sim_singlesp_genome_${INDEX} --single ~/scratch/mt_sim/ld_7point5mill_sample_processed.fasta -l 200 -s 5 -t 30 > log_singlesp_genome_sim_${INDEX}.txt 2>&1
	&& mv sim_singlesp_genome_${INDEX}/abundance.tsv ../sim_singlesp_genome_${INDEX}_abundance.txt
done
echo "second real genome"
indexes=("2_1" "2_2" "2_3" "2_4" "2_5" "2_6" "2_7")
FASTQ="SRR769395"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	cd ~/scratch/ensemble_bact/part_${INDEX}/
	time kallisto index -i singlesp_genome_${INDEX}_index *.dna.genome.fa > log_index_${INDEX}.txt 2>&1 \
	&& time kallisto quant -i singlesp_genome_${INDEX}_index -o real_singlesp_genome_${INDEX} ~/scratch/mt_real/${FASTA}_1.fastq.gz ~/scratch/mt_real/${FASTA}_2.fastq.gz -t 30 > log_align_${FASTQ}_${INDEX}.txt 2>&1 \
	&& mv kallisto_singlesp_genome_${INDEX}/abundance.tsv ../kallisto_singlesp_genome_${INDEX}_abundance.txt
done
