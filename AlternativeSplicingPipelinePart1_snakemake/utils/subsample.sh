#!/bin/bash

#SBATCH --job-name=TKU4354_subsampling
#SBATCH --mem=120G
#SBATCH --partition=pe2
#SBATCH --output=./logs/TKU4354_subsampling.log
#SBATCH --error=./logs/TKU4354_subsampling.err

module load seqtk
cd /gpfs/commons/groups/landau_lab/rraviram/Suva_lab_GBM/Splicing_ONT/ONT_Splicing_TKU4354/input_files/1.ONT_fastq/

seqtk sample TKU4354.fastq 5000000 > TKU4354_subsample_5b.fastq
