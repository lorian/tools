#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT path info or file extension (which needs to be .mfa)
DIR=~/scratch/clark_$FASTA/
mkdir $DIR
cd $DIR
mkdir Custom
cd Custom \
&& python ~/tools/M_prep_clark_files_single.py ~/scratch/${FASTA}.mfa \
&& cd ~/scratch/clark \
&& time ./set_targets.sh $DIR custom \
&& time ./classify_metagenome.sh -P ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq -R clark_${FASTA} \
&& time ./estimate_abundance.sh -F clark_${FASTA}.csv -D $DIR > ${FASTA}.clark \
