#!/bin/bash
indexes=("entero" "klebsiella" "shigella" "misc" "ecoli_1" "ecoli_2" "ecoli_3" "ecoli_4")
FASTQ="SRR769395"
for INDEX in "${indexes[@]}"
do
	echo -e "\n\e[7m$INDEX\e[0m"
	cd ~/scratch/ensembl_cdna/real_${INDEX}/
	time kallisto index -i real_${INDEX}_index *.fa > ../log_index_${INDEX}.txt 2>&1 \
	&& time kallisto quant -i real_${INDEX}_index -o real_${INDEX} ~/scratch/mt_real/${FASTQ}_1.fastq.gz ~/scratch/mt_real/${FASTQ}_2.fastq.gz -t 30 > ../log_align_${FASTQ}_${INDEX}.txt 2>&1 \
	&& mv real_${INDEX}/abundance.tsv ../real_${INDEX}_abundance.txt
done
