## bam_check.R script 
# runs samtools quickcheck to confirm the listed BAM files
# exist and/or are readable by the current environment and confirms the BAM
# files have an intact file header this will not catch incorrect genome builds,
# read coverage, or incomplete files/corrupted files
args = commandArgs(trailingOnly=TRUE)
bam <- snakemake@input[["bam"]]
outname <- snakemake@output[[1]]

source("scripts/funcs.R")
bamCheck(x = bam)