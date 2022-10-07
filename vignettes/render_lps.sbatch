#!/bin/bash

#SBATCH --job-name=render
#SBATCH --output=render_%A_%a.out
#SBATCH --error=render_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --array=1-112

module load R


CHROMS=('chr1' 'chr1' 'chr1' 'chr1' 'chr1' 'chr1' 'chr1' 'chr1' 'chr1'
 'chr10' 'chr10' 'chr10' 'chr10' 'chr10' 'chr10'
 'chr11' 'chr11' 'chr11' 'chr11' 'chr11' 'chr11' 'chr11'
 'chr12' 'chr12' 'chr12' 'chr12' 'chr12' 'chr12' 
 'chr13' 'chr13' 
 'chr14' 'chr14'
 'chr15' 'chr15' 'chr15' 'chr15' 'chr15'
 'chr16' 'chr16' 'chr16'
 'chr17' 'chr17' 'chr17' 'chr17' 'chr17' 'chr17' 'chr17'
 'chr18' 'chr18'
 'chr19' 'chr19' 'chr19' 'chr19'
 'chr2' 'chr2' 'chr2' 'chr2' 'chr2' 'chr2' 'chr2' 'chr2' 'chr2'
 'chr20' 'chr20' 'chr20' 'chr20'
 'chr21' 'chr21'
 'chr22' 'chr22'
 'chr3' 'chr3' 'chr3' 'chr3' 'chr3' 'chr3' 'chr3' 
 'chr4' 'chr4' 'chr4' 'chr4'
 'chr5' 'chr5' 'chr5' 'chr5' 'chr5' 'chr5'
 'chr6'	'chr6' 'chr6' 'chr6' 'chr6' 'chr6' 'chr6' 'chr6'
 'chr7' 'chr7' 'chr7' 'chr7' 'chr7'
 'chr8' 'chr8' 'chr8' 'chr8' 'chr8' 'chr8'
 'chr9' 'chr9' 'chr9' 'chr9' 'chr9'
 'chrX' 'chrX'
)

CHUNKS=(1 2 3 4 5 6 7 8 9
 1 2 3 4 5 6
 1 2 3 4 5 6 7
 1 2 3 4 5 6
 1 2
 1 2
 1 2 3 4 5
 1 2 3
 1 2 3 4 5 6 7
 1 2
 1 2 3 4
 1 2 3 4 5 6 7 8 9
 1 2 3 4
 1 2
 1 2
 1 2 3 4 5 6 7
 1 2 3 4
 1 2 3 4 5 6
 1 2 3 4 5 6 7 8
 1 2 3 4 5
 1 2 3 4 5 6
 1 2 3 4 5
 1 2
)

CHROM=${CHROMS[${SLURM_ARRAY_TASK_ID}]}
CHUNK=${CHUNKS[${SLURM_ARRAY_TASK_ID}]}

Rscript -e "rmarkdown::render('loop_strength_lps_loops.Rmd', params = list(chrom = '$CHROM', chunk = $CHUNK))"
