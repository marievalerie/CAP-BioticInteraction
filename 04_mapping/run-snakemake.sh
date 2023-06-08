#! /bin/bash
#
#$ -S /bin/bash
#$ -N snakejob
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load anaconda3/2020.02
conda activate snakemake
module load bowtie2/2.3.5.1

a=1; b=$NSLOTS
THREADS=$((b-a)) #snakemake requires 1 thread for job monitoring

snakemake --cores ${THREADS} --snakefile Snakefile_mapping --use-conda --rerun-incomplete
