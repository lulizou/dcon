#!/bin/bash

#SBATCH --job-name=renderx
#SBATCH --output=renderx_%A_%a.out
#SBATCH --error=renderx_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --array=0-7

module load R/4.1

STARTS=('0' '10e6' '20e6' '30e6' '40e6' '50e6' '60e6' '70e6' '80e6' '90e6' '100e6'
	'110e6' '120e6' '130e6' '140e6' '150e6' '160e6')
ENDS=('10e6' '20e6' '30e6' '40e6' '50e6' '60e6' '70e6' '80e6' '90e6' '100e6'
       	'110e6' '120e6' '130e6' '140e6' '150e6' '160e6' '170e6')
DFS=(50 100 250 500 750 1000 2500 5000)
START=${STARTS[${SLURM_ARRAY_TASK_ID}]}
END=${ENDS[${SLURM_ARRAY_TASK_ID}]}
DF=${DFS[${SLURM_ARRAY_TASK_ID}]}

Rscript -e "rmarkdown::render('rdsprite-deconvolve-xist.Rmd', params = list(start = '60e6', end='150e6', df=$DF, mol='L1Md_DNA_matrix_count_RNA_10000_clusters_2-100000_chrchrX.mtx'))"
