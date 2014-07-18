# Version 1.1 -- used for testM

#sed -e '/>/s/^/@/' -e '/>/s/$/#/' Martin_etal_TextS3_13Dec2011_poundless.fasta | tr -d "\n" | tr "@" "\n" | sort -t "|" -k1n | tr "#" "\n" | sed -e '/^$/d' > Martin_etal_TextS3_13Dec2011_sorted.fasta \
#&& python M_process_metagenome_reference.py Martin_etal_TextS3_13Dec2011_sorted.fasta \
#&&bowtie2-build Martin_etal_TextS3_13Dec2011_sorted_genomes_0.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_0 \
#&& bowtie2-build Martin_etal_TextS3_13Dec2011_sorted_genomes_1.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_1 \
#&& bowtie2-build Martin_etal_TextS3_13Dec2011_sorted_genomes_2.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_2 \
#&& bowtie2-build Martin_etal_TextS3_13Dec2011_sorted_plasmids_0.mfa Martin_etal_TextS3_13Dec2011_sorted_plasmids_0 \

bowtie2 -a -t -p 30 --local -x Martin_etal_TextS3_13Dec2011_sorted_genomes_0 -1 illumina_100species.1.fq -2 illumina_100species.2.fq -S testM_g0.SAM \
&& bowtie2 -a -t -p 30 --local -x Martin_etal_TextS3_13Dec2011_sorted_genomes_1 -1 illumina_100species.1.fq -2 illumina_100species.2.fq -S testM_g1.SAM \
&& bowtie2 -a -t -p 30 --local -x Martin_etal_TextS3_13Dec2011_sorted_genomes_2 -1 illumina_100species.1.fq -2 illumina_100species.2.fq -S testM_g2.SAM \
&& bowtie2 -a -t -p 30 --local -x Martin_etal_TextS3_13Dec2011_sorted_plasmids_0 -1 illumina_100species.1.fq -2 illumina_100species.2.fq -S testM_p.SAM \
&& python M_merge_sam_files.py testM_g0.SAM testM_g1.SAM testM_g2.SAM testM_p.SAM testE.SAM \
&& mv combined_file.sam testE_M.sam \
&& samtools view -bhS testE_M.sam > testE_M.bam \
&& ulimit -n 5000 \
&& samtools sort -n testE_M.bam testE_M_sorted \

#&& cat illumina_100genomes.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_0.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_1.mfa Martin_etal_TextS3_13Dec2011_sorted_genomes_2.mfa Martin_etal_TextS3_13Dec2011_sorted_plasmids_0.mfa > Martin_etal_TextS3_13Dec2011_sorted_i100.mfa \

&& express -f 0.85 --max-indel-size 100 -B 20 Martin_etal_TextS3_13Dec2011_sorted_i100.mfa testE_M_sorted.bam
