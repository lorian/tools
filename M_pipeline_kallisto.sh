#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT file extension (which needs to be .mfa)
time /usr/local/bin/kallisto_fast_loose index -i {$FASTA}_index {$FASTA}.mfa \
&& time /usr/local/bin/kallisto_fast_loose quant --fast-and-loose -i {$FASTA}_index -o $FASTA illumina_100species_trimmed.1.fq illumina_100species_trimmed.2.fq \
&& cd $FASTA \
&& mv abundance.tsv ../{$FASTA}_abundance.txt
