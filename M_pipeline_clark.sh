#!/bin/bash
FASTA=$1 #accepts argument of fasta name, WITHOUT path info or file extension (which needs to be .mfa)
DIR=~/scratch/clark_$FASTA/
mkdir $DIR
cd $DIR
mkdir Custom
cd Custom \
&& echo python ~/tools/M_prep_clark_files_single.py ~/scratch/${FASTA}.mfa \
&& python ~/tools/M_prep_clark_files_single.py ~/scratch/${FASTA}.mfa \
&& cd ~/scratch/clark \
&& echo ./set_targets.sh $DIR custom \
&& time ./set_targets.sh $DIR custom \
&& echo ./classify_metagenome.sh -P ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq -R clark_${FASTA} \
&& time ./classify_metagenome.sh -P ~/scratch/kallisto_data/illumina_100species_trimmed.1.fq ~/scratch/kallisto_data/illumina_100species_trimmed.2.fq -R clark_${FASTA} \
&& echo ./estimate_abundance.sh -F clark_${FASTA}.csv -D $DIR > ${FASTA}.clark \
&& time ./estimate_abundance.sh -F clark_${FASTA}.csv -D $DIR > ${FASTA}.clark \
&& echo Output in ${FASTA}.clark
