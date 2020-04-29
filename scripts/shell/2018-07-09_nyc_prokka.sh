#!/bin/bash
#SBATCH -n 8
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -o 2018-07-09-nyc-prokka.out
#SBATCH -e 2018-07-09-nyc-prokka.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"
#SBATCH --array=0-950

module load perl/5.24.0
module load java/jdk-1.8u112

export PATH=$PATH:/n/data1/hsph/immid/grad/software/prokka-1.13/bin/

readarray -t f < good_samples.txt

prokka --force --outdir ${f[${SLURM_ARRAY_TASK_ID}]} --locustag ${f[${SLURM_ARRAY_TASK_ID}]} --prefix ${f[${SLURM_ARRAY_TASK_ID}]} --genus Neisseria --species gonorrhoeae --strain ${f[${SLURM_ARRAY_TASK_ID}]} --cpus 8 ../spades/assemblies/${f[${SLURM_ARRAY_TASK_ID}]}_contigs_filtered.fa
