#!/bin/bash
indexes=(1 2 3 4 5)
fasta = "SRR769395"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	time kallisto quant -i singlesp_genome_${INDEX}_index -o real_singlesp_genome_${INDEX} ~/scratch/mt_real/${FASTA}_1.fastq.gz ~/scratch/mt_real/${FASTA}_2.fastq.gz -t 30 > log_align_${FASTA}_${INDEX}.txt 2>&1 \
	&& mv real_singlesp_genome_${INDEX}/abundance.tsv ../real_singlesp_genome_${FASTA}_${INDEX}_abundance.txt
done
