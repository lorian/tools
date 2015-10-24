#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT path info or file extension (which needs to be .mfa)
DIR=~/scratch/kraken_$FASTA
./kraken-build --download-taxonomy --db $DIR \
&& python ~/tools/M_prep_kraken_files.py ~/scratch/${FASTA}.mfa \
&& ./kraken-build --add-to-library ~/scratch/${FASTA}_kraken.fasta --db $DIR \
&& time ./kraken-build --build --db $DIR \
&& time ./kraken --db $DIR --fastq-input --paired ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq >$
&& time ./kraken-translate --db $DIR ${FASTA}.kraken > labeled_${FASTA}.kraken \
