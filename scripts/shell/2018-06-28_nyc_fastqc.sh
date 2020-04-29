#!/bin/bash
#SBATCH -n 1
#SBATCH --mem-per-cpu=1G
#SBATCH	-p short
#SBATCH -t 20:00
#SBATCH -o 2018-06-28-nyc-fastqc.out
#SBATCH -o 2018-06-28-nyc-fastqc.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-2763

module load java/jdk-1.8u112

readarray f < fastqs.txt

mkdir -p 2018-06-28-nyc-fastqc/

/home/tm241/grad/software/FastQC-0.11.7/FastQC/fastqc -q -t 1 -o 2018-06-28-nyc-fastqc/ -d /home/tm241/tm241/tmp/ ${f[${SLURM_ARRAY_TASK_ID}]}
