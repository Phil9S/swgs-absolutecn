args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(dplyr))

bam_in <- as.character(args[1])
sample_name <- as.character(args[2])
meta <- as.character(args[3])
rds <- as.character(args[4])
outdir <- args[5]
bin <- as.numeric(args[6])
project <- args[7]

fit.qc <- read.table(file = meta,header = T,sep = "\t",na.strings = "")

relative_smoothed <- readRDS(rds)
read.data <- phenoData(relative_smoothed)@data

fit.qc.filt <- fit.qc %>%
  filter(use == TRUE)

if(!sample_name %in% fit.qc.filt$SAMPLE_ID){ 
  paste0("Sample ids: ",sample_name," had no fits passing QC")
  quit(save="no")
}

fit.qc.filt$total.reads <- read.data$total.reads[match(x = fit.qc.filt$SAMPLE_ID,read.data$name)]
fit.qc.filt$ratio <- round(fit.qc.filt$downsample_depth / fit.qc.filt$total.reads,digits = 2)

perc <- fit.qc.filt %>%
   filter(SAMPLE_ID == sample_name) %>%
   .$ratio
  
if( perc <= 0.96){
  cmd.downsample <- paste("samtools view -s ", perc," -b ",bam_in," > ",paste0(outdir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/downsampled_bams/",sample_name,".bam"))
  cmd.index <- paste0("samtools index ",paste0(outdir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/downsampled_bams/",sample_name,".bam"))
   
  system(cmd.downsample)
  system(cmd.index)
    
 }else{
  cmd.copy <- paste0("cp ",bam_in," ",paste0(outdir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/downsampled_bams/",sample_name,".bam"))
  system(cmd.copy)
 }
