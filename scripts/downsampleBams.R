args = commandArgs(trailingOnly=TRUE)

bam_in <- snakemake@input[["bam"]]
meta <- snakemake@input[["meta"]]
rds <- snakemake@input[["rds"]]
outdir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
outname <- snakemake@output[[1]]
sample_name <- snakemake@params[["sample"]]
prplpu <- as.logical(snakemake@params[["prplpu"]])
filetype <- snakemake@params[["filetype"]]
reference <- snakemake@params[["reference"]]

options(scipen=999)

fit.qc <- read.table(file = meta,header = T,sep = "\t",na.strings = "")
fit.qc.filt <- fit.qc %>%
  dplyr::filter(SAMPLE_ID == sample_name) %>%
  dplyr::filter(use == TRUE)

if(prplpu){
  # Get read count using SAMtools if prplpu is set for all samples as the 
  # stage_1 fitting is not performed
  cmd.totalreads <- paste0("samtools view -c -F 260 ",bam_in)
  tot.reads <- as.numeric(system(cmd.totalreads,intern=TRUE))
  read.data <- data.frame(name = fit.qc.filt$SAMPLE_ID,total.reads = tot.reads)
} else {
  relative_smoothed <- readRDS(rds)
  read.data <- Biobase::phenoData(relative_smoothed)@data
}

fit.qc.filt$total.reads <- read.data$total.reads[match(x = fit.qc.filt$SAMPLE_ID,read.data$name)]
fit.qc.filt$ratio <- round(fit.qc.filt$downsample_depth / fit.qc.filt$total.reads,digits = 4)

perc <- fit.qc.filt %>% .$ratio

# Notes on downsample ratio
# - Read ratio is greater than 0.96 (i.e close to original or higher than 
# available reads) are symlink of original.
# - CRAM files will always be downsampled into bams regardless of ratio between 
# downsample depth to prevent repeat decompression to BAM
# - Total reads Error catch for read ratios lower than 1e-4 (Large bin / high 
# coverages / Low ploidy / High purity) implemented to prevent too few reads 
# being sampled.

if(perc > 1){
  perc <- 1
} else if(perc == 0){
  message("downsample ratio too low - set to 0.0001 percent")
  perc <- 0.0001
}

if(filetype == "CRAM"){
  cmd.downsample <- paste("samtools view -s ", perc," -T ",reference," -b ",
                          bam_in," > ",outname)
  cmd.index <- paste0("samtools index ",outname)

  system(cmd.downsample)
  system(cmd.index)
 
} else {
  if( perc <= 0.96){
    cmd.downsample <- paste("samtools view -s ", perc," -b ",
                            bam_in," > ",outname)
    cmd.index <- paste0("samtools index ",outname)
   
    system(cmd.downsample)
    system(cmd.index)
    
  }else{
    cmd.copy <- paste0("ln -s ",bam_in," ",outname)
    cmd.index <- paste0("samtools index ",outname)
    system(cmd.copy)
    system(cmd.index)
  }
}
