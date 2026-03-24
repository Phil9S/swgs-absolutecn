## Functions and vars used by swgs-abscn

# define vars
`%>%` <- dplyr::`%>%`

fittingColumnNames <- c("SAMPLE_ID","ploidy","purity","clonality","rmse",
                        "downsample_depth","powered","TP53cn","expected_TP53_AF",
                        "homozygousLoss","MedianSegVar")

dsFittingColumnNames <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","precPloidy",
                          "precPurity","TP53freq","clonality.pre","clonality.post", 
                          "rmse.pre","rmse.post","paired.ends","total.reads",
                          "used.reads","expected.variance","loess.span","loess.family",
                          "downsample_depth","powered","TP53cn.pre","TP53cn.post",
                          "expected_TP53_AF.pre","expected_TP53_AF.post","smooth",
                          "homozygousLoss","MedianSegVar.pre","MedianSegVar.post",
                          "rank_clonality","pl_diff","new_state_n","new_state","use","notes")


# Define helper functions
# Performs quickcheck of BAM files
bamCheck <- function(x=NULL,outname=NULL){
  if(is.null(x)){
    stop("no BAMs provided")
  }
  
  log_vector <- c()
  check_vector <- c()
  for(i in x){
    ## Call samtools using cmdline
    cmd <- paste0("samtools quickcheck -q ",i)
    if(system(cmd) == 0){
      log_vector <- append(log_vector, paste0("BAM valid - ",i))
      check_vector <- append(check_vector,FALSE)
    } else {
      log_vector <- append(log_vector,paste0("BAM invalid or missing - ",i))
      check_vector <- append(check_vector,TRUE)
    }
  }
  
  if(any(check_vector)){
    outname <- gsub(pattern = "ok",replacement = "invalid",outname)
    writeLines(text = as.character(log_vector),con = outname)	
  } else {
    writeLines(text = as.character(log_vector),con = outname)
  }
}

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

calculateSegmentVar <- function(abs_seg=NULL,abs_cn=NULL){
  segRLE <- rle(abs_seg)
  
  segVar <- c()
  for(i in 1:length(segRLE$lengths)){
    if(i == 1){
      strt_idx <- 1
      end_idx <- segRLE$lengths[i]
      segVar <- append(segVar,var(abs_cn[strt_idx:end_idx]))
    } else {
      start_idx <- max(cumsum(segRLE$lengths[1:i-1]))
      strt_idx <- start_idx + 1
      if(i == length(segRLE$lengths)){
        end_idx <- segRLE$lengths[i] + strt_idx - 1
      } else {
        end_idx <- segRLE$lengths[i] + strt_idx
      }
      segVar <- append(segVar,var(abs_cn[strt_idx:end_idx]))
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

getSegTable <- function(x){
  if(inherits(x,what = "QDNAseqCopyNumbers",which = F)){
    sn <- Biobase::assayDataElement(x,"segmented")
    fd <- Biobase::fData(x)
    use <- fd$use
    fdfiltfull <- fd[use,]
    sn <- sn[use,]
    if(is.null(ncol(sn))){
      sampleName <- Biobase::sampleNames(x)
      sn <- as.data.frame(sn)
      colnames(sn) <- sampleName
    }
    segTable <- c()
    for(s in colnames(sn)){
      for(c in unique(fdfiltfull$chromosome)){
        snfilt <- sn[fdfiltfull$chromosome == c,colnames(sn) == s]
        fdfilt <- fdfiltfull[fdfiltfull$chromosome == c,]
        sn.rle <- rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        
        segtmp <- do.call(rbind,lapply(1:length(sn.rle$lengths), function(s){
          from <- fdfilt$start[starts[s]]
          to <- fdfilt$end[ends[s]]
          segValue <- sn.rle$value[s]
          c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }))
        
        segTableRaw <- cbind(segtmp,sample = rep(s,times=nrow(segtmp)))
        segTable <- rbind(segTable,segTableRaw)
      }
    }
    segTable <- as.data.frame(segTable)
    colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
    segTable$segVal <- round(as.numeric(segTable$segVal),3)
    segTable$segVal[segTable$segVal < 0] <- 0
    segTable$start <- as.numeric(segTable$start)
    segTable$end <- as.numeric(segTable$end)
    return(segTable)
  } else {
    # NON QDNASEQ BINNED DATA FUNCTION
    stop("segtable error")
  }
}

plotProfile <- function(relcn=NULL,ploidy=NA,purity=NA,clonality=NA,rmse=NA,yrange=NULL){
  if(is.null(relcn)){
    stop("no data")
  }
  # Y axis range
  if(is.null(yrange)){
    if(ploidy > 5){
      yrange = 15
    } else {
      yrange = 10
    }
  } else{
    if(!is.numeric(yrange)){
      stop("yrange must be a integer")
    }
  }
  # Plot abs fit
  sub <- paste0("ploidy=",round(ploidy,2)," | ",
                " purity=",round(purity,2)," | ",
                " MAE=",round(clonality,3)," | ",
                " RMSE=",round(rmse,3))
  
  QDNAseqmod::plot(relcn,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,
                   ylim=c(0,yrange),ylab="Absolute tumour CN")
  abline(h=1:yrange-1, col = "blue")
  mtext(sub,side = 1,line = 3.5)
}

filterFitTable <- function(table = NULL,metadata = NULL,filter_underpowered = NULL,filter_homozygous = NULL,af_cutoff = NULL,ranks = 10){
  
  fitTable <- dplyr::left_join(table,metadata,by="SAMPLE_ID") %>%
    dplyr::select(-file) %>%
    dplyr::relocate(PATIENT_ID,.after = SAMPLE_ID) %>%
    dplyr::relocate(TP53freq,smooth,.after = expected_TP53_AF) %>%
    dplyr::relocate(precPloidy,precPurity,.after = purity)
  
  ## Apply hard filters
  ##  filter under powered fits when config variable is TRUE
  if(filter_underpowered){
    fitTable <- fitTable %>%
      dplyr::filter(powered == 1) 
  }
  ## filter high prop homozygous loss when config variable is TRUE
  if(filter_homozygous){
    fitTable <- fitTable %>%
      dplyr::filter(homozygousLoss <= homozygous_prop)
  }
  
  # standard filtering
  filtered_results <- fitTable %>%
    dplyr::group_by(SAMPLE_ID, ploidy) %>%
    dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% #rank clonality within a unique ploidy state 
    dplyr::filter(rank_clonality == 1) %>% #select ploidy with the lowest clonality within a unique ploidy state 
    dplyr::group_by(SAMPLE_ID) %>%
    dplyr::top_n(-ranks, wt = clonality) %>% # select top 10 ploidy states with the lowest clonality values
    dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% # rank by clonality within a sample across ploidy in top 10
    # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
    dplyr::filter(is.na(TP53freq) | dplyr::near(expected_TP53_AF,TP53freq,tol = af_cutoff)) %>%  
    dplyr::arrange(PATIENT_ID, SAMPLE_ID)
  
  pruned_results <- filtered_results %>%
    dplyr::arrange(SAMPLE_ID, ploidy) %>%
    dplyr::group_by(SAMPLE_ID) %>%
    dplyr::mutate(pl_diff = abs(ploidy - dplyr::lag(ploidy))) %>% #, pu_diff = abs(purity - dplyr::lag(purity)) not used
    dplyr::mutate(new_state_n = dplyr::row_number() == 1 | pl_diff > 0.3) %>%
    dplyr::mutate(new_state = cumsum(new_state_n)) %>%
    dplyr::group_by(SAMPLE_ID, new_state) %>%
    dplyr::filter(rank_clonality == min(rank_clonality)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(use = rep(NA,times=nrow(.)),notes = rep(NA,times=nrow(.))) %>%
    dplyr::arrange(ploidy)
  
  return(list(filtered = filtered_results,pruned = pruned_results))
}
