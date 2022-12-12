#!/bin/bash

#SBATCH --job-name=renderx
#SBATCH --output=renderx_%A_%a.out
#SBATCH --error=renderx_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --array=0-21

module load R/4.1

DF=250
START="60e6"
END="120e6"


shopt -s nullglob
array=(/rafalab/lzou/rdsprite/dna_matrix/*mtx.mtx)
f=${array[${SLURM_ARRAY_TASK_ID}]}

Rscript -e "rmarkdown::render('rdsprite-deconvolve-xist.Rmd', params = list(start = '60e6', end='120e6', df=$DF, mol='$f'))"
