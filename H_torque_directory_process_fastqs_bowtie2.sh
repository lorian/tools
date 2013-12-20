#PBS -35 nodes=1:ppn=1,mem=100g,walltime=168:00:00
#PBS -q batch
#PBS -V

# qsub -v BASENAME="DH3_hESN-termGFP_line1",SUFFIX="_a" H_torque_directory_process_fastqs.sh
# $BASENAME is the basename of the fastq file
# $SUFFIX is the suffix to be used to distinguish runs

cd /home/lorian/hockemeyer

../bowtie2-2.1.0/bowtie2 -t -p35 --rdg6,5 --rfg6,5 --score-minL,-.6,-.4 -xfused_ensembl_updated_cdna -U${BASENAME}.fastq -S${BASENAME}${SUFFIX}.SAM

