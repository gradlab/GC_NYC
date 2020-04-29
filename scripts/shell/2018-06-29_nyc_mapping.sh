#!/bin/bash
#SBATCH -n 8
#SBATCH --mem-per-cpu=1G
#SBATCH	-p short
#SBATCH -t 10:00
#SBATCH -o 2018-06-29-nyc-mapping.out
#SBATCH -e 2018-06-29-nyc-mapping.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-1381

module load gcc/6.2.0
module load java/jdk-1.8u112
module load bwa/0.7.15
module load samtools/1.3.1
module load picard/2.8.0

readarray -t f < samples.txt

bwa mem -M -t 8 /n/data1/hsph/immid/grad/Tatum/references/fastas/Neisseria_gonorrhoeae_NCCP11945.fa "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_1.fastq.gz" "/n/data1/hsph/immid/grad/gonococcus/NYC_data/fastq_files/${f[${SLURM_ARRAY_TASK_ID}]}_2.fastq.gz" | samtools sort -o "${f[${SLURM_ARRAY_TASK_ID}]}.bam" -
java -Xmx4g -jar $PICARD/picard-2.8.0.jar MarkDuplicates I="${f[${SLURM_ARRAY_TASK_ID}]}.bam" O="${f[${SLURM_ARRAY_TASK_ID}]}.marked.bam" M="${f[${SLURM_ARRAY_TASK_ID}]}.metrics"
