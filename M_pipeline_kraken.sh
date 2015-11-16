#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT path info or file extension (which needs to be .mfa)
DIR=~/scratch/kraken_$FASTA
echo "./kraken-build --download-taxonomy --db $DIR" \
&& ./kraken-build --download-taxonomy --db $DIR \
&& echo "python ~/tools/M_prep_kraken_files.py ~/scratch/${FASTA}.mfa" \
&& python ~/tools/M_prep_kraken_files.py ~/scratch/${FASTA}.mfa \
&& echo "./kraken-build --add-to-library ~/scratch/${FASTA}_kraken.fasta --db $DIR" \
&& ./kraken-build --add-to-library ~/scratch/${FASTA}_kraken.fasta --db $DIR \
&& echo "./kraken-build --build --db $DIR" \
&& time ./kraken-build --build --db $DIR \
&& echo "./kraken --db $DIR --fastq-input --paired ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq > ${FASTA}.kraken" \
&& time ./kraken --db $DIR --fastq-input --paired ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq > ${FASTA}.kraken \
&& echo "./kraken-translate --db $DIR ${FASTA}.kraken > labeled_${FASTA}.kraken" \
&& time ./kraken-translate --db $DIR ${FASTA}.kraken > labeled_${FASTA}.kraken \
&& echo "Output in labeled_${FASTA}.kraken"
