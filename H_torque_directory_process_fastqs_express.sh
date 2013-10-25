#PBS -l nodes=1:ppn=1,mem=2g,walltime=168:00:00
#PBS -q batch
#PBS -V

# qsub -v BASENAME="DH3_hESN-termGFP_line1",SUFFIX="_a" H_torque_directory_process_fastqs.sh
# $BASENAME is the basename of the sam file
# $SUFFIX is the suffix to be used to distinguish runs

cd /home/lorian/hockemeyer

express -f0.75 -B20 -o${BASENAME}${SUFFIX} fused_ensembl_updated_cdna.fa ${BASENAME}${SUFFIX}.SAM
