#!/bin/bash
indexes=("2_1" "2_2" "2_3" "2_4" "2_5" "2_6")
FASTQ="SRR769395"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	cd ~/scratch/ensembl_bact/part_${INDEX}/
	time kallisto quant -i singlesp_genome_${INDEX}_index -o real_singlesp_genome_${INDEX} ~/scratch/mt_real/${FASTQ}_1.fastq.gz ~/scratch/mt_real/${FASTQ}_2.fastq.gz -t 30 > log_align_${FASTQ}_${INDEX}.txt 2>&1 \
	&& mv kallisto_singlesp_genome_${INDEX}/abundance.tsv ../kallisto_singlesp_genome_${INDEX}_abundance.txt
done
