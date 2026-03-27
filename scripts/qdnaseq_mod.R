## qdnaseq_mod.R
# Run QDNASeqmod to generate segmented ReadCount CN profiles for
# use downstream This script performs both pre and post downsampling profile
# binning and segmentation, includes parameters for PE/SE data, and
# auto-smoothing functions
args = commandArgs(trailingOnly=TRUE)
bin.size <- as.numeric(snakemake@params[["bin"]])
pairedEnd <- as.logical(snakemake@params[["pairedEnd"]])
sample <- snakemake@params[["sample"]]
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
use_seed <- snakemake@params[["use_seed"]]
seed_val <- snakemake@params[["seed_val"]]
bam_list <- snakemake@input[["bams"]]
outname <- snakemake@output[["rds"]]
genome <- as.character(snakemake@params[["genome"]])
autoSmooth <- as.logical(snakemake@params[["autoSmooth"]])
smoothThreshold <- as.numeric(snakemake@params[["smoothThreshold"]])
# Samples to smooth
smooth <- unique(metadata$smooth[metadata$SAMPLE_ID == sample])
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
# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)
## generate annotation file either by preloading calculated files or generating new one
bins <- QDNAseqmod::getBinAnnotations(binSize=bin.size,genome=genome)
readCounts <- QDNAseqmod::binReadCounts(bamfiles = bam_list,bins = bins,
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
copyNumbersSegmentedSmooth <- smooth_sample(relcn = copyNumbersSegmented,
                                            smooth=smooth,maxSegs=maxSegs,
                                            seed=seed)

# save output to file
saveRDS(copyNumbersSegmentedSmooth,outname)
