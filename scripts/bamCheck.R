## bam_check.R
# runs samtools quickcheck to confirm the listed BAM files
# exist and/or are readable by the current environment and confirms the BAM
# files have an intact file header this will not catch incorrect genome builds,
# read coverage, or incomplete files/corrupted files
args = commandArgs(trailingOnly=TRUE)

## Set args based on if using snakemake or cmdline
if(exists("snakemake")){
  bam <- snakemake@input[["bam"]]
  fileType <- snakemake@params[["filetype"]]
  okOut <- snakemake@output[["ok"]]
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "bamCheck"))
  
  bam <- opts$bam
  fileType <- opts$fileType
  okOut <- opts$okOut
}

rswgsabsolutecn::bamCheck(x = bam,filetype = fileType,outname = okOut)
