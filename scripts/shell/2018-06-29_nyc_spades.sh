#!/bin/bash
#SBATCH -n 8
#SBATCH --mem-per-cpu=4G
#SBATCH	-p short
#SBATCH -t 3:00:00
#SBATCH -o 2018-06-29-nyc-spades.out
#SBATCH -e 2018-06-29-nyc-spades.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=2-1381


readarray -t f < samples.txt

export PATH=$PATH:/n/data1/hsph/immid/grad/software/SPAdes-3.12.0-Linux/bin/


spades.py -t 8 --careful -1  "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_1.fastq.gz" -2 "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_2.fastq.gz" -o "${f[${SLURM_ARRAY_TASK_ID}]}"
spades.py -t 8 --careful --plasmid --pe1-1  "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_1.fastq.gz" --pe1-2 "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_2.fastq.gz" -o "${f[${SLURM_ARRAY_TASK_ID}]}_plasmid"
