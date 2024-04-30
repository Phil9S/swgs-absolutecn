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
prplpu <- snakemake@params[["prplpu"]]

#print(prplpu)

fit.qc <- read.table(file = meta,header = T,sep = "\t",na.strings = "")
fit.qc.filt <- fit.qc %>%
  filter(SAMPLE_ID == sample_name) %>%
  filter(use == TRUE)

if(prplpu == "TRUE"){
  cmd.totalreads <- paste0("samtools view -c -F 260 ",bam_in)
  tot.reads <- as.numeric(system(cmd.totalreads,intern=TRUE))
  read.data <- data.frame(name=fit.qc.filt$SAMPLE_ID,total.reads=tot.reads)
  print(read.data)
} else {
  relative_smoothed <- readRDS(rds)
  read.data <- phenoData(relative_smoothed)@data
}

fit.qc.filt$total.reads <- read.data$total.reads[match(x = fit.qc.filt$SAMPLE_ID,read.data$name)]
fit.qc.filt$ratio <- round(fit.qc.filt$downsample_depth / fit.qc.filt$total.reads,digits = 3)

perc <- fit.qc.filt %>%
   .$ratio

# If read ratio is greater than 0.96 (i.e close to original or higher than available reads)
if( perc <= 0.96){
  cmd.downsample <- paste("samtools view -s ", perc," -b ",bam_in," > ",outname)
  cmd.index <- paste0("samtools index ",outname)
   
  system(cmd.downsample)
  system(cmd.index)
    
 }else{
  cmd.copy <- paste0("cp ",bam_in," ",outname)
  cmd.index <- paste0("samtools index ",outname)
  system(cmd.copy)
  system(cmd.index)
 }
