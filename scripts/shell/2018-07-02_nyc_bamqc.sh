#!/bin/bash
#SBATCH -n 8
#SBATCH --mem-per-cpu=1G
#SBATCH	-p short
#SBATCH -t 10:00
#SBATCH -o 2018-07-02-nyc-bamqc.out
#SBATCH -e 2018-07-02-nyc-bamqc.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-1381

module load java/jdk-1.8u112

readarray -t f < samples.txt

~/grad/software/qualimap_v2.2.1/qualimap bamqc -nt 8 -bam "../mapping/${f[${SLURM_ARRAY_TASK_ID}]}.marked.bam" -outdir "${f[${SLURM_ARRAY_TASK_ID}]}"
