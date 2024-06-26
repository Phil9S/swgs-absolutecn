#Rscript qdnaseq.R -m metadata.csv" -b 30 -n 5
args = commandArgs(trailingOnly=TRUE)

bin.size <- as.numeric(snakemake@params[["bin"]])
#ncores <- as.numeric(snakemake@resources[["cpus"]])
ncores <- 1
output_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
metafile <- snakemake@input[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bam_list <- snakemake@input[["bam"]]
outname <- snakemake@output[[1]]
use_seed <- snakemake@params[["use_seed"]]
seed_val <- snakemake@params[["seed_val"]]


suppressMessages(library(parallel))
suppressMessages(library(tidyverse))
suppressMessages(library(Biobase))
suppressMessages(library(QDNAseqmod))
suppressMessages(library(plyr))

metadata <- metadata[metadata$use == "TRUE",]

sampleIds <- unique(metadata$SAMPLE_ID)

## generate annotation file either by preloading calculated files or generating new one
bins <- getBinAnnotations(binSize=bin.size)

# Samples to smooth
smoothed_samples <- as.character(metadata$SAMPLE_ID[metadata$smooth == "TRUE"])

readCounts <- mclapply(X=bam_list, FUN=binReadCounts, bins=bins, mc.cores=ncores,chunkSize=1e7)
## if copyNumbersSegment file exists read it else generate it
# apply filter based on loess fit residuals and encode/1000-genome balcklist
readCountsFiltered <- mclapply(X=readCounts, FUN=applyFilters, mc.cores=1)

# estimate correction for GC content and mappability
readCountsFiltered <- mclapply(X=readCountsFiltered, FUN=estimateCorrection, mc.cores=1)
# apply the correction for GC content and mappability
copyNumbers <- mclapply(X=readCountsFiltered, FUN=correctBins, mc.cores=1)

#bring back to readcount space 
assayDataElement(copyNumbers[[1]],"copynumber") <- assayDataElement(copyNumbers[[1]],"copynumber") * median(assayDataElement(readCountsFiltered[[1]], "fit"), na.rm=T)

# smooth outliers (Data is now ready to be analyzed with a downstream package of choice (exportBins))
copyNumbersSmooth <- mclapply(X=copyNumbers, FUN=smoothOutlierBins, mc.cores=1)

# Implement seeding to prevent variable segments on repeated runs
if(use_seed){
  seed <- as.character(seed_val)
} else {
  seed <- NULL
}

# perform segmentation on bins and save it
copyNumbersSegmented <- mclapply(X=copyNumbersSmooth, FUN=segmentBins, transformFun="sqrt", mc.cores=ncores,seeds=seed)

changeSampleName <- function(CNsObj)
  {
    sampleNames(CNsObj) <- gsub(x=unlist(lapply(strsplit(pData(CNsObj)$name,split="\\."),function(x) x[1])),
					pattern="_ds",
					replacement="")
    pData(CNsObj)$name <- sampleNames(CNsObj)
    return(CNsObj)
  }

#copyNumbersSegmented <- mclapply(X=copyNumbersSegmented, FUN=changeSampleName, mc.cores=ncores)

# For each
smooth_samples <- function(obj){
  # Index and subselect sample
  #ind <- which(colnames(copyNumbersSegmentedSmooth)==sample)
  relcn <- obj
  # Check if smoothing needed
  smooth.bool <- FALSE
  relative_tmp <- NULL
  if(sampleNames(obj) %in% smoothed_samples){
    smooth.bool <- TRUE
    currsamp <- relcn
    maxseg<-300
    sdadjust<-1.5
    condition <- fData(currsamp)$use
    segments <- assayDataElement(currsamp, "segmented")[condition, , drop=FALSE]
    segments<-apply(segments,2,rle)
    segnum<-as.numeric(lapply(segments,function(x){length(x$lengths)}))
    while(sum(segnum>maxseg)&sdadjust<5)
    {
      currsamp<-segmentBins(currsamp, transformFun="sqrt",undo.SD=sdadjust,seeds=seed)
      segments <- assayDataElement(currsamp, "segmented")[condition, , drop=FALSE]
      segments<-apply(segments,2,rle)
      segnum<-as.numeric(lapply(segments,function(x){length(x$lengths)}))
      sdadjust<-sdadjust+0.5
    }
    relative_tmp <- currsamp
    relcn <- relative_tmp
  }
  return(relcn)
}

copyNumbersSegmentedSmooth <- mclapply(X=copyNumbersSegmented, FUN=smooth_samples, mc.cores=ncores)

if(is.na(pData(object=copyNumbersSegmentedSmooth[[1]])$loess.span)){
        stop(paste0(sampleNames(copyNumbersSegmented)," BAM failed loess fitting. Remove this file from sample sheet"))
}


saveRDS(copyNumbersSegmentedSmooth,outname)
