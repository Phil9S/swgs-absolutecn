## downsampleBams.R
# Notes on downsample ratio
# - Read ratio is greater than 0.96 (i.e close to original or higher than 
# available reads) are symlink of original.
# - CRAM files will always be downsampled into bams regardless of ratio between 
# downsample depth to prevent repeat decompression to BAM
# - Total reads Error catch for read ratios lower than 1e-4 (Large bin / high 
# coverages / Low ploidy / High purity) implemented to prevent too few reads 
# being sampled.
args = commandArgs(trailingOnly=TRUE)

bam <- snakemake@input[["bam"]]
meta <- snakemake@input[["meta"]]
rds <- snakemake@input[["rds"]]
bamOut <- snakemake@output[["bam"]]
sampleName <- snakemake@params[["sample"]]
prplpu <- as.logical(snakemake@params[["prplpu"]])
fileType <- snakemake@params[["filetype"]]
reference <- snakemake@params[["reference"]]

options(scipen=999)
`%>%` <- dplyr::`%>%`

fitQC <- read.table(file = meta,header = T,sep = "\t",na.strings = "")
fitQCFilt <- fitQC %>%
  dplyr::filter(SAMPLE_ID == sampleName) %>%
  dplyr::filter(use == TRUE)

if(prplpu){
  # Get read count using SAMtools if prplpu is set for all samples as the 
  # stage_1 fitting is not performed
  cmdTotalreads <- paste0("samtools view -c -F 260 ",bam)
  totReads <- as.numeric(system(cmdTotalreads,intern=TRUE))
  readData <- data.frame(name = fitQCFilt$SAMPLE_ID,total.reads = totReads)
} else {
  relative_smoothed <- readRDS(rds)
  readData <- Biobase::phenoData(relative_smoothed)@data
}

fitQCFilt$total.reads <- readData$total.reads[match(x = fitQCFilt$SAMPLE_ID,readData$name)]
fitQCFilt$ratio <- round(fitQCFilt$downsample_depth / fitQCFilt$total.reads,digits = 4)

perc <- fitQCFilt %>% .$ratio

if(perc > 1){
  perc <- 1
} else if(perc == 0){
  message("downsample ratio too low - set to 0.0001 percent")
  perc <- 0.0001
}

if(fileType == "CRAM"){
  if(perc <= 0.96){
    cmdDownsample <- paste("samtools view -s ",perc," -T ",reference," -b ",bam," > ",bamOut)
  } else {
    cmdDownsample <- paste("samtools view -T ",reference," -b ",bam," > ",bamOut)
  }
  cmdIndex <- paste0("samtools index ",bamOut)

  system(cmdDownsample)
  system(cmdIndex)
 
} else {
  if(perc <= 0.96){
    cmdDownsample <- paste("samtools view -s ", perc," -b ",bam," > ",bamOut)
    cmdIndex <- paste0("samtools index ",bamOut)
   
    system(cmdDownsample)
    system(cmdIndex)
    
  } else {
    cmdCopy <- paste0("ln -s ",bam," ",bamOut)
    cmdIndex <- paste0("samtools index ",bamOut)
    system(cmdCopy)
    system(cmdIndex)
  }
}

