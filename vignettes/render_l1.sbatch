#!/bin/bash

#SBATCH --job-name=render
#SBATCH --output=render_%A_%a.out
#SBATCH --error=render_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --array=0-19

module load R/4.1

CHROMS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 
        'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chrX')
CHROM=${CHROMS[${SLURM_ARRAY_TASK_ID}]}

Rscript -e "rmarkdown::render('l1_observed_expected.Rmd', params = list(chrom = '$CHROM'))"
