python ~/tools/M_mash_dist_parallel.py kraken_SRS_trimmed.fastq \
&& cat mash_real_kraken_SRS_trimmed_* > mash_real_kraken_trimmed.txt \
&& python ~/tools/M_mash_kallisto_pipeline.py ~/scratch/realdata/mash_real_kraken_trimmed.txt 50 ~/scratch/realdata/mash_SRS/ \
&& time kallisto index -i real_trimmed_index ~/scratch/realdata/mash_SRS/*.cat.fa \
&& time kallisto quant -i real_trimmed_index ~/scratch/realdata/kraken_trimmed_all.1.fastq ~/scratch/realdata/kraken_trimmed_all.2.fastq -o kallisto_real_trimmed -t 15 \
&& mv kallisto_real_trimmed/abundance.tsv kallisto_real_trimmed_abundance.txt
