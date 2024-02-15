#!/bin/bash

#SBATCH --mem=32g
#SBATCH --job-name=submit_annotation

module load R/3.6.0

workingdir=$1
#threedata=$1
#fivedata=$2
#output=$3
run_files=$2

Rscript "$run_files"/DTU_junction_annotation.new.R $workingdir

