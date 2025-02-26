#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=32g 
#SBATCH --job-name=combine_output

run_files=$1
outputsdir=$2
cell_meta_data=$3
finalout=$4

Rscript "$run_files"/merge_final_output_comb_patient_merge_counts.R $outputsdir $cell_meta_data $finalout


