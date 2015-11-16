#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT path info or file extension (which needs to be .mfa)
echo kallisto index -i ${FASTA}_index ~/scratch/${FASTA}.mfa \
&& time kallisto index -i ${FASTA}_index ~/scratch/${FASTA}.mfa \
&& echo kallisto quant -i ${FASTA}_index -o $FASTA illumina_100species_trimmed.1.fq illumina_100species_trimmed.2.fq \
&& time kallisto quant -i ${FASTA}_index -o $FASTA illumina_100species_trimmed.1.fq illumina_100species_trimmed.2.fq \
&& cd $FASTA \
&& mv abundance.tsv ../${FASTA}_abundance.txt \
&& cd ../ \
&& echo Output in ${FASTA}_abundance.txt

