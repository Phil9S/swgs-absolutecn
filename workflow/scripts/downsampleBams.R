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

# Redirect R output, messages, warnings and errors to the Snakemake rule log
if(exists("snakemake")){
  logFile <- file(snakemake@log[[1]], open = "wt")
  sink(logFile, type = "output")
  sink(logFile, type = "message")
}


if(exists("snakemake")){
  bam <- snakemake@input[["bam"]]
  meta <- snakemake@input[["meta"]]
  rds <- snakemake@input[["rds"]]
  bamOut <- snakemake@output[["bam"]]
  sampleName <- snakemake@params[["sample"]]
  prplpu <- as.logical(snakemake@params[["prplpu"]])
  fileType <- snakemake@params[["filetype"]]
  reference <- snakemake@params[["reference"]]
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "downsampleBams"))
  
  bam <- opts$bam
  meta <- opts$meta
  rds <- opts$rds
  bamOut <- opts$bamOut
  sampleName <- opts$sampleName
  prplpu <- opts$prplpu
  fileType <- opts$fileType
  reference <- opts$reference
}

options(scipen=999)
`%>%` <- dplyr::`%>%`

fitQC <- read.table(file = meta,header = T,sep = "\t",na.strings = "",
  colClasses = c(SAMPLE_ID = "character", PATIENT_ID = "character"))
fitQCFilt <- fitQC %>%
  dplyr::filter(SAMPLE_ID == sampleName) %>%
  dplyr::filter(use == "TRUE")

# Exactly one selected fit is required per sample to set the downsample depth
if(nrow(fitQCFilt) != 1){
  stop(paste0(sampleName," - expected 1 fit with use == TRUE in ",meta,
              " but found ",nrow(fitQCFilt)))
}

# Shell-quoted paths for use in system() commands (paths may contain spaces)
bamQ <- shQuote(bam)
bamOutQ <- shQuote(bamOut)
referenceQ <- shQuote(reference)

if(prplpu){
  # Get read count using SAMtools if prplpu is set for all samples as the
  # stage_1 fitting is not performed
  if(fileType == "CRAM"){
    cmdTotalreads <- paste0("samtools view -T ",referenceQ," -c -F 260 ",bamQ)
  } else {
    cmdTotalreads <- paste0("samtools view -c -F 260 ",bamQ)
  }
  totReads <- suppressWarnings(as.numeric(system(cmdTotalreads,intern=TRUE)))
  if(length(totReads) != 1 || is.na(totReads)){
    stop(paste0("Failed to count reads: ", cmdTotalreads))
  }
  readData <- data.frame(name = fitQCFilt$SAMPLE_ID,total.reads = totReads)
} else {
  relative_smoothed <- readRDS(rds)
  readData <- Biobase::phenoData(relative_smoothed)@data
}

fitQCFilt$total.reads <- readData$total.reads[match(x = fitQCFilt$SAMPLE_ID,readData$name)]
fitQCFilt$ratio <- round(fitQCFilt$downsample_depth / fitQCFilt$total.reads,digits = 4)

perc <- fitQCFilt %>% .$ratio

if(is.na(perc)){
  stop(paste0(sampleName," - downsample ratio is NA (downsample_depth: ",
              fitQCFilt$downsample_depth,", total.reads: ",fitQCFilt$total.reads,")"))
}

if(perc > 1){
  perc <- 1
} else if(perc == 0){
  message("downsample ratio too low - set to 0.01 percent")
  perc <- 0.0001
}

if(fileType == "CRAM"){
  if(perc <= 0.96){
    cmdDownsample <- paste0("samtools view --subsample ",perc," --subsample-seed 99 -T ",referenceQ," -b -o ",bamOutQ," ",bamQ)
  } else {
    cmdDownsample <- paste0("samtools view -T ",referenceQ," -b -o ",bamOutQ," ",bamQ)
  }
} else {
  if(perc <= 0.96){
    cmdDownsample <- paste0("samtools view --subsample ",perc," --subsample-seed 99 -b -o ",bamOutQ," ",bamQ)
  } else {
    cmdDownsample <- paste0("cp -P ",bamQ," ",bamOutQ)
  }
}
cmdIndex <- paste0("samtools index ",bamOutQ)

# Stop on failure so a truncated BAM is not passed
status <- system(cmdDownsample)
if(status != 0){
  unlink(bamOut)
  stop(paste0("Downsampling failed (exit status ",status,"): ",cmdDownsample))
}

status <- system(cmdIndex)
if(status != 0){
  unlink(c(bamOut,paste0(bamOut,".bai")))
  stop(paste0("BAM indexing failed (exit status ",status,"): ",cmdIndex))
}

