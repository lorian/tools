#!/bin/bash
indexes=(2 3 4 5)
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	cd part_${INDEX}/
	time kallisto index -i singlesp_genome_${INDEX}_index *.dna.genome.fa \
	&& time kallisto quant -i singlesp_genome_${INDEX}_index -o kallisto_singlesp_genome_${INDEX} --single ~/scratch/mt_sim/ld_7point5mill_sample_processed.fasta -l 200 -s 5 -t 30 \
	&& mv kallisto_singlesp_genome_${INDEX}/abundance.tsv ../kallisto_singlesp_genome_${INDEX}_abundance.txt
done
