#!/bin/bash

#SBATCH --job-name=test
#SBATCH --mem=10G
#SBATCH --partition=pe2
#SBATCH --output=test.log
#SBATCH --error=test.err
#SBATCH --cpus-per-task=5

module load anaconda3/10.19
module load java/1.9
module load samtools
module load bedtools
module load racon

source activate use_mamba
source activate sicelore2.0 

snakemake --profile profile_snakemake
