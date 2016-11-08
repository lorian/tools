#!/bin/bash
indexes=(1 2 3 4 6 7)
FASTQ="SRR769516"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	echo ${FASTQ}
	time kallisto quant -i singlesp_genome_${INDEX}_index -o real_DNA_singlesp_genome_${INDEX} ~/scratch/mt_real/${FASTA}_1.fastq ~/scratch/mt_real/${FASTA}_2.fastq -t 30 \
	&& mv real_DNA_singlesp_genome_${INDEX}/abundance.tsv ../real_DNA_singlesp_genome_${FASTQ}_${INDEX}_abundance.txt
done
