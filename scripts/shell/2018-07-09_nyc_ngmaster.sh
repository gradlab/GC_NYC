#!/bin/bash
#SBATCH -n 1
#SBATCH --mem-per-cpu=4G
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -o 2018-06-01-ngmaster-NYC.out
#SBATCH -e 2018-06-01-ngmaster-NYC.err
#SBATCH --mail-type=END
#SBATCH --mail-user=mortimer@hsph.harvard.edu

export PATH=$PATH:/n/data1/hsph/immid/grad/software/ngmaster/isPcr_x86_64/:/n/data1/hsph/immid/grad/software/ngmaster/
ngmaster --printseq 2018-07-09_nyc_ngmast_alleles.fa ../spades/assemblies/*filtered.fa > 2018-07-09_nyc_ngmast.txt
