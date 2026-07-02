## bam_check.R
# runs samtools quickcheck to confirm the listed BAM files
# exist and/or are readable by the current environment and confirms the BAM
# files have an intact file header this will not catch incorrect genome builds,
# read coverage, or incomplete files/corrupted files
args = commandArgs(trailingOnly=TRUE)
bam <- snakemake@input[["bam"]]
filetype <- snakemake@params[["filetype"]]
outname <- snakemake@output[["ok"]]

source(file.path(snakemake@scriptdir,"funcs.R"))
bamCheck(x = bam,filetype = filetype,outname = outname)
