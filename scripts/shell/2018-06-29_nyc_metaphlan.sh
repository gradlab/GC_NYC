#!/bin/bash
#SBATCH -n 2
#SBATCH --mem-per-cpu=2G
#SBATCH	-p short
#SBATCH -t 30:00
#SBATCH -o 2018-06-29-nyc-metaphlan.out
#SBATCH -e 2018-06-29-nyc-metaphlan.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-1381

module load gcc/6.2.0
module load bowtie2/2.2.9

readarray -t f < samples.txt
mkdir -p 2018-06-29-nyc-metaphlan/

zcat "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_1.fastq.gz" "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_2.fastq.gz" | /n/data1/hsph/immid/grad/software/metaphlan2/metaphlan2.py --bowtie2out 2018-06-29-nyc-metaphlan/${f[${SLURM_ARRAY_TASK_ID}]}.bowtie2.bz2 --nproc 2 --input_type fastq -o 2018-06-29-nyc-metaphlan/metaphlan2_${f[${SLURM_ARRAY_TASK_ID}]}.txt

rm 2018-06-29-nyc-metaphlan/${f[${SLURM_ARRAY_TASK_ID}]}.bowtie2.bz2
