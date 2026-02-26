## Functions used by swgs-abscn
#define helper functions
# converts readdepth to copy number given purity and single copy depth
depthtocn<-function(x,purity,seqdepth){
  (x/seqdepth-2*(1-purity))/purity
}

get_gene_seg <- function(target=NULL,abs_data=NULL){
  to_use <- Biobase::fData(abs_data)$use
  cn_obj <- abs_data[to_use,]
  bin_pos <- Biobase::fData(cn_obj)[,c("chromosome","start","end")]
  position <- as.numeric(stringr::str_split(string = target,pattern = ":|-",simplify = T))
  mapply(assign,c("chr","start","end"),position,MoreArgs=list(envir=parent.frame()))
  
  gene_chr_pos <- bin_pos[bin_pos$chromosome == chr,]
  min_start <- min(which(min(abs(gene_chr_pos$start - start)) == abs(gene_chr_pos$start - start)))
  min_end <- max(which(min(abs(gene_chr_pos$end - end)) == abs(gene_chr_pos$end - end)))
  if(gene_chr_pos$start[min_start] > start & min_start != 1){
    min_start <- min_start - 1
  }
  if(gene_chr_pos$end[min_end] < end & min_end != length(gene_chr_pos$end)){
    min_end <- min_end + 1
  }
  index_min <- which(bin_pos$chromosome == chr & bin_pos$start == gene_chr_pos[min_start,2])
  index_max <- which(bin_pos$chromosome == chr & bin_pos$end == gene_chr_pos[min_end,3])
  gene_pos <- seq.int(index_min,index_max,1)
  return(gene_pos)
}

getDownsampleDepth <- function(ploidy=NULL,purity=NULL,nbins_ref_genome=NULL,rpc=15,ratio=1.098901){
  # original implementation
  #(((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)
  if(any(is.null(ploidy),is.null(purity),is.null(nbins_ref_genome))){
    stop("missing parameters")
  }
  
  cellploidy <- purity * ploidy + (2 * (1-purity))
  relratio <- (cellploidy/(ploidy*purity)) / purity
  
  readRatio <- relratio * rpc * cellploidy * nbins_ref_genome * ratio
  return(readRatio)
}

calculateSegmentVar <- function(abs_cn=NULL,abs_cnbin=NULL){
  segRLE <- rle(abs_cn)
  
  segVar <- c()
  for(i in 1:length(segRLE$lengths)){
    if(i == 1){
      strt_idx <- 1
      end_idx <- segRLE$lengths[i]
      segVar <- append(segVar,var(abs_cnbin[strt_idx:end_idx]))
    } else {
      start_idx <- max(cumsum(segRLE$lengths[1:i-1]))
      strt_idx <- start_idx + 1
      if(i == length(segRLE$lengths)){
        end_idx <- segRLE$lengths[i] + strt_idx - 1
      } else {
        end_idx <- segRLE$lengths[i] + strt_idx
      }
      segVar <- append(segVar,var(abs_cnbin[strt_idx:end_idx]))
    }
  }
  
  medianVar <- median(segVar)
  return(medianVar)
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

collapse_rds <- function(rds.list){
  comb <- rds.list[[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]]
      comb <- Biobase::combine(comb,add)
    }
    rds.obj <- comb
  } else {
    rds.obj <- comb
  } 
  return(rds.obj)
}