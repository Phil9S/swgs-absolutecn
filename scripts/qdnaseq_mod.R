#Rscript Run QDNASeqmod to generate segmented ReadCount CN profiles
args = commandArgs(trailingOnly=TRUE)

## Load snakemake args
bin.size <- as.numeric(snakemake@params[["bin"]])
output_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
sample <- snakemake@params[["sample"]]
metafile <- snakemake@params[["meta"]]
use_seed <- snakemake@params[["use_seed"]]
seed_val <- snakemake@params[["seed_val"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bam_list <- snakemake@input[["bams"]]
outname <- snakemake@output[[1]]
genome <- as.character(snakemake@params[["genome"]])

# Samples to smooth
#smoothed_samples <- as.character(metadata$SAMPLE_ID[metadata$smooth == "TRUE"])
smooth <- metadata$smooth[metadata$SAMPLE_ID == sample]

# Implement seeding to prevent variable segments on repeated runs
if(use_seed){
  seed <- as.character(seed_val)
} else {
  seed <- NULL
}

# segment smoothing function
# Adjusted recursively to drop max segments below threshold using maximum StdDev 
# difference in means during segmentation splits in CBS
smooth_sample <- function(relcn=NULL,smooth=FALSE,maxSegs=300,seed=NULL){
  
  if(is.null(relcn)){
    stop("segment smoothing provided with no data")
  }
  
  stopifnot(is.logical(smooth))
  stopifnot(is.numeric(maxSegs),maxSegs > 22)
  
  # Check if smoothing needed
  relative_tmp <- NULL
  if(smooth){
    currsamp <- relcn
    
    maxseg <- maxSegs
    sdadjust <- 1.5
    
    condition <- Biobase::fData(currsamp)$use
    
    segments <- Biobase::assayDataElement(currsamp, "segmented")[condition, , drop=FALSE]
    segments <- apply(segments,2,rle)
    segnum<-as.numeric(lapply(segments,function(x){length(x$lengths)}))
    
    while(segnum > maxseg & sdadjust < 5){
      currsamp <- QDNAseqmod::segmentBins(currsamp, transformFun="sqrt",undo.SD=sdadjust,seeds=seed)
      
      segments <- Biobase::assayDataElement(currsamp, "segmented")[condition, , drop=FALSE]
      segments <- apply(segments,2,rle)
      segnum <- as.numeric(lapply(segments,function(x){length(x$lengths)}))
      
      sdadjust <- sdadjust + 0.5
    }
    #relative_tmp <- currsamp
    relcn <- currsamp
  }
  return(relcn)
}

## generate annotation file either by preloading calculated files or generating new one
bins <- QDNAseqmod::getBinAnnotations(binSize=bin.size,genome=genome)

#readCounts <- mclapply(X=bam_list, FUN=binReadCounts, bins=bins,mc.cores=ncores,chunkSize=1e7)
##### THIS IS MAKING SOME EXTERNAL CONNECTION ON FIRST USE
readCounts <- QDNAseqmod::binReadCounts(bamfiles = bam_list,bins = bins,chunkSize = 1e7)
# apply filter based on loess fit residuals and encode/1000-genome blacklist
readCountsFiltered <- QDNAseqmod::applyFilters(object = readCounts)
# estimate correction for GC content and mappability
readCountsFiltered <- QDNAseqmod::estimateCorrection(object = readCountsFiltered)

## Edge case error check for BAM files which pass a quick check but are insufficent to 
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
copyNumbersSegmented <- QDNAseqmod::segmentBins(object = copyNumbersSmooth,transformFun="sqrt",seeds=seed)

# smooth copy number segmentation
copyNumbersSegmentedSmooth <- smooth_sample(relcn = copyNumbersSegmented,smooth=smooth,maxSegs=300,seed=seed)

# save output to file
saveRDS(copyNumbersSegmentedSmooth,outname)
