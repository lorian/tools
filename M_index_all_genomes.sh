#!/bin/bash
indexes=(1 2 3 4 5 6 7)
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	cd ../part_${INDEX}/
	time kallisto index -i singlesp_genome_${INDEX}_index *.dna.genome.fa > log_index_${INDEX}.txt 2>&1 \
	&& time kallisto quant -i singlesp_genome_${INDEX}_index -o kallisto_singlesp_genome_${INDEX} --single ~/scratch/mt_sim/ld_7point5mill_sample_processed.fasta -l 200 -s 5 -t 30 > log_align_${INDEX}.txt 2>&1 \
	&& mv kallisto_singlesp_genome_${INDEX}/abundance.tsv ../kallisto_singlesp_genome_${INDEX}_abundance.txt
done
