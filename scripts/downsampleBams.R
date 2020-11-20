args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))

bam_in <- snakemake@input[["bam"]]
meta <- snakemake@input[["meta"]]
rds <- snakemake@input[["rds"]]
outdir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
outname <- snakemake@output[[1]]
sample_name <- snakemake@params[["sample"]]

fit.qc <- read.table(file = meta,header = T,sep = "\t",na.strings = "")

relative_smoothed <- readRDS(rds)
read.data <- phenoData(relative_smoothed)@data

fit.qc.filt <- fit.qc %>%
  filter(use == TRUE)

fit.qc.filt$total.reads <- read.data$total.reads[match(x = fit.qc.filt$SAMPLE_ID,read.data$name)]
fit.qc.filt$ratio <- round(fit.qc.filt$downsample_depth / fit.qc.filt$total.reads,digits = 2)

perc <- fit.qc.filt %>%
   filter(SAMPLE_ID == sample_name) %>%
   .$ratio

if( perc <= 0.96){
  cmd.downsample <- paste("samtools view -s ", perc," -b ",bam_in," > ",outname)
  cmd.index <- paste0("samtools index ",outname)
   
  system(cmd.downsample)
  system(cmd.index)
    
 }else{
  cmd.copy <- paste0("cp ",bam_in," ",outname)
  system(cmd.copy)
 }
