#!/bin/bash
#SBATCH -n 1
#SBATCH --mem-per-cpu=50M
#SBATCH	-p short
#SBATCH -t 1:00
#SBATCH -o 2018-06-29-nyc-pseudogenomes.out
#SBATCH -e 2018-06-29-nyc-pseudogenomes.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-950

readarray -t f < good_samples.txt

./pilonVCFtoFasta.py ../variants/"${f[${SLURM_ARRAY_TASK_ID}]}_pilon.vcf"
