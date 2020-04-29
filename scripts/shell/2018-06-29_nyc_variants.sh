#!/bin/bash
#SBATCH -n 1
#SBATCH --mem-per-cpu=16G
#SBATCH	-p short
#SBATCH -t 30:00
#SBATCH -o 2018-06-29-nyc-mapping.out
#SBATCH -e 2018-06-29-nyc-mapping.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-1381

module load java/jdk-1.8u112
module load samtools/1.3.1

readarray -t f < samples.txt

samtools index ../mapping/"${f[${SLURM_ARRAY_TASK_ID}]}.marked.bam" ../mapping/"${f[${SLURM_ARRAY_TASK_ID}]}.marked.bam.bai" 

java -Xmx16g -jar /n/data1/hsph/immid/grad/software/pilon/pilon-1.16.jar --genome /n/data1/hsph/immid/grad/Tatum/references/fastas/Neisseria_gonorrhoeae_NCCP11945.fa --frags ../mapping/"${f[${SLURM_ARRAY_TASK_ID}]}.marked.bam" --output "${f[${SLURM_ARRAY_TASK_ID}]}_pilon" --variant --vcf --mindepth 10 --minmq 20
