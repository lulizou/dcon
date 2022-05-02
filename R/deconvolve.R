library(argparse)
library(mgcv)
library(dplyr)
library(ggplot2)
library(reshape2)
library(splines)
library(survival)
library(spatstat)
library(data.table)
set.seed(1337)
options(scipen=999)
source('/n/irizarryfs01/lzou/dcon/R/utils.R')


parser <- ArgumentParser()
parser$add_argument('-i', '--input', type='character', required=T,
                    help='Input counts in .bedgraph format')
parser$add_argument('-c', '--chrom', type='character', default='chr22',
                    help='Which chromosome to do')
parser$add_argument('-w', '--window', type='integer', default=1e5,
                    help='Size of window to perform deconvolution in')
parser$add_argument('-o', '--output', type='character', required=T,
                    help='output file prefix')
args <- parser$parse_args()


hichip <- fread(args$input, header=F, 
                col.names=c('chrom','start','end','count')) %>% 
  filter(chrom==params$chrom)
cuts <- fread('/n/irizarryfs01/lzou/resources/hg38_mboi.bed', header=F, 
              col.names=c('chrom','start','end','name','idk1','strand')) %>% 
  filter(chrom==params$chrom)
cuts$round <- floor(cuts$start/500)*500
cut_var <- cuts %>% group_by(round) %>% summarise(cut=n())
hic <- readRDS('/n/irizarryfs01/lzou/salva/K562_HiC_',args$chrom,'_500.rds')
hichip$pos <- seq(1,nrow(hichip))






