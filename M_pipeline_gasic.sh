#!/bin/bash
BOOTSTRAP=$1 #accepts argument of # of boostraps
echo python ~/tools/M_directory_create_bowtie2_indexes.py ncbi_data \
&& python ~/tools/M_directory_create_bowtie2_indexes.py ncbi_data \
&& echo python gasic-r16/run_mappers_nonduplicate.py ncbi_data_names.txt illumina_100species.all.fq -m bowtie2 -i ncbi_data/%s -o i100_gasic/%s.SAM \
&& time python gasic-r16/run_mappers_nonduplicate.py ncbi_data_names.txt illumina_100species.all.fq -m bowtie2 -i ncbi_data/%s -o i100_gasic/%s.SAM \
&& export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3.4/dist-packages \
&& echo python gasic-r16/create_matrix.py -s mason_illumina -r ncbi_data/%s.mfa -m bowtie2 -i ncbi_data/%s -t ~/metagenome/tmp/ -o i100_gasic/i100_matrix.npy ncbi_data_names.txt \
&& time python gasic-r16/create_matrix.py -s mason_illumina -r ncbi_data/%s.mfa -m bowtie2 -i ncbi_data/%s -t ~/metagenome/tmp/ -o i100_gasic/i100_matrix.npy ncbi_data_names.txt \
&& echo python gasic-r16/correct_abundances.py -m i100_gasic/i100_matrix.npy -s i100_gasic/%s.SAM -b $BOOTSTRAP -o gasic_output_${BOOTSTRAP}.txt ncbi_data_names.txt ]
&& time python gasic-r16/correct_abundances.py -m i100_gasic/i100_matrix.npy -s i100_gasic/%s.SAM -b $BOOTSTRAP -o gasic_output_${BOOTSTRAP}.txt ncbi_data_names.txt \
&& echo Output in gasic_output_${BOOTSTRAP}.txt
