#!/bin/bash
#SBATCH -n 12
#SBATCH --mem-per-cpu=8G
#SBATCH -p priority
#SBATCH -t 8:00:00
#SBATCH -o 2019-01-23-nyc-gubbins.out
#SBATCH -e 2019-01-23-nyc-gubbins.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu
#SBATCH --constraint="scratch2"

run_gubbins.py --threads 12 --prefix 2019-01-23_nyc_metadataComplete ../2019-01-23_nyc_pseudogenomes_metadataComplete.fasta
