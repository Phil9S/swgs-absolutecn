## qdnaseq_mod.R
# Run QDNASeqmod to generate segmented ReadCount CN profiles for
# use downstream This script performs both pre and post downsampling profile
# binning and segmentation, includes parameters for PE/SE data, and
# auto-smoothing functions
args = commandArgs(trailingOnly=TRUE)

if(exists("snakemake")){
  
  bin <- as.numeric(snakemake@params[["bin"]])
  sampleName <- snakemake@params[["sample"]]
  meta <- snakemake@params[["meta"]]
  bam <- snakemake@input[["bam"]]
  rdsOut <- snakemake@output[["rds"]]
  genome <- as.character(snakemake@params[["genome"]])
  autoSmooth <- as.logical(snakemake@params[["autoSmooth"]])
  smoothThreshold <- as.numeric(snakemake@params[["smoothThreshold"]])
  use_seed <- snakemake@params[["use_seed"]]
  seed_val <- snakemake@params[["seed_val"]]
  pairedEnd <- as.logical(snakemake@params[["pairedEnd"]])

} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "qdnaseqMod"))
  
  bin <- as.numeric(opts$bin)
  sampleName <- opts$sampleName
  meta <- opts$meta
  bam <- opts$bam
  rdsOut <- opts$rdsOut
  genome <- as.character(opts$genome)
  autoSmooth <- as.logical(opts$autoSmooth)
  smoothThreshold <- as.numeric(opts$smoothThreshold)
  use_seed <- opts$
  seed_val <- opts$
  pairedEnd <- as.logical(opts$pairedEnd)  
}

options(scipen = 999)
metadata <- read.table(file = meta,header=T,sep="\t")

# Samples to smooth
smooth <- unique(metadata$smooth[metadata$SAMPLE_ID == sampleName])

# Allow for automatic smoothing above threshold
if(autoSmooth){
  smooth <- autoSmooth
  maxSegs <- smoothThreshold
} else {
  maxSegs <- 300
}
# Implement seeding to prevent variable segments on repeated runs
if(use_seed){
  seed <- as.character(seed_val)
} else {
  seed <- NULL
}

## generate annotation file either by preloading calculated files or generating new one
bins <- QDNAseqmod::getBinAnnotations(binSize=bin,genome=genome)
readCounts <- QDNAseqmod::binReadCounts(bamfiles = bam,bins = bins,
                                        chunkSize = 1e7,pairedEnds = pairedEnd)

# apply filter based on loess fit residuals and encode/1000-genome blacklist
readCountsFiltered <- QDNAseqmod::applyFilters(object = readCounts)

# estimate correction for GC content and mappability
readCountsFiltered <- QDNAseqmod::estimateCorrection(object = readCountsFiltered)

## Edge case error check for BAM files which pass a quick check but are insufficient to 
# produce a loess model from the data present. No way to fix this other than removing the file.
if(is.na(Biobase::pData(readCountsFiltered)$loess.span)){
  stop(paste0(Biobase::sampleNames(readCountsFiltered),
              " BAM failed loess fitting. Remove this file from analysis"))
}

# apply the correction for GC content and mappability
copyNumbers <- QDNAseqmod::correctBins(object = readCountsFiltered)

# bring back to readcount space
medianRC <- median(Biobase::assayDataElement(readCountsFiltered, "fit"), na.rm=T)
Biobase::assayDataElement(copyNumbers,"copynumber") <- Biobase::assayDataElement(copyNumbers,"copynumber") * medianRC

# smooth outlier bins
copyNumbersSmooth <- QDNAseqmod::smoothOutlierBins(object = copyNumbers)

# perform segmentation on bins and save it
copyNumbersSegmented <- QDNAseqmod::segmentBins(object = copyNumbersSmooth,
                                                transformFun="sqrt",seeds=seed)

# smooth copy number segmentation 
copyNumbersSegmentedSmooth <- rswgsabsolutecn::smoothSample(relcn = copyNumbersSegmented,
                                            smooth=smooth,maxSegs=maxSegs,
                                            seed=seed)

# save output to file
saveRDS(copyNumbersSegmentedSmooth,rdsOut)
