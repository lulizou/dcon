#!/bin/bash

#SBATCH --job-name=k562
#SBATCH --output=k562_%A_%a.out
#SBATCH --error=k562_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=10G
#SBATCH --array=0-21

module load R/4.1

CHROMS=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 
        'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')
CHROM=${CHROMS[${SLURM_ARRAY_TASK_ID}]}

Rscript -e "rmarkdown::render('loop_strength_k562_mnase_ctcf.Rmd', params = list(chrom = '$CHROM'))"
