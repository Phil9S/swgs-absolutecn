## Added by PS
args = commandArgs(trailingOnly=TRUE)
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
filter_underpowered <- as.logical(snakemake@params[["filter_underpowered"]]) # Filter for powered fits only
filter_homozygous <- as.logical(snakemake@params[["filter_homozygous"]]) # filter homozygous loss
homozygous_prop <- as.numeric(snakemake@params[["homozygous_prop"]])

#threads
cores <- as.numeric(snakemake@threads)
#cores <- 8
doMC::registerDoMC(cores)

# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)
outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_PRE_down_sampling/")

# read in relative CN data
# collapse rds files function
rds.filename <- snakemake@input[["rds"]]
rds.list <- lapply(rds.filename,FUN=function(x){readRDS(x)})
# Combine and load rds objects - COULD CHANGE TO one-by-one rather than jointly
relative_smoothed <- collapse_rds(rds.list)
saveRDS(relative_smoothed,paste0(outpath,project,"_",bin,"kb_relSmoothedCN.rds"))

filelist <- snakemake@input[["cl"]]
fitTable <- do.call(rbind,
			lapply(filelist,FUN = function(x){
			  ## ADD SAMPLE NAME IN PREVIOUS STEP RATHER THAN FROM FILE NAME
				#n <- gsub(pattern="_clonality.tsv",rep="",x=x)
        #prefix <- paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",project,"_")
        #n <- gsub(pattern=prefix,rep="",x=n)
				tab <- read.table(x,sep="\t",skip=1)
				#tab <- cbind(rep(n,times=nrow(tab)),tab)
				return(tab)
			}))
colnames(fitTable) <- fittingColumnNames

`%>%` <- dplyr::`%>%`
fitTable <- dplyr::left_join(fitTable,metadata,by="SAMPLE_ID") %>%
  dplyr::select(-file) %>%
  dplyr::relocate(PATIENT_ID,.after = SAMPLE_ID) %>%
  dplyr::relocate(TP53freq,smooth,.after = expected_TP53_AF) %>%
  dplyr::relocate(precPloidy,precPurity,.after = purity)
  #select(SAMPLE_ID,PATIENT_ID,ploidy,purity,clonality,downsample_depth,powered,TP53cn,expected_TP53_AF,TP53freq,smooth,homozygousLoss)

## Apply hard filters
##  filter underpowered fits when config variable is TRUE
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
  #dplyr::select(SAMPLE_ID, PATIENT_ID, everything()) %>% redundant select
  dplyr::group_by(SAMPLE_ID, ploidy) %>%
  dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% #rank clonality within a unique ploidy state 
  dplyr::filter(rank_clonality == 1) %>% #select ploidy with the lowest clonality within a unique ploidy state 
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::top_n(-10, wt = clonality) %>% # select top 10 ploidy states with the lowest clonality values
  dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% # rank by clonality within a sample across ploidy in top 10
  # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
  dplyr::filter(is.na(TP53freq) | dplyr::near(expected_TP53_AF,TP53freq,tol = af_cutoff)) %>%  
  dplyr::arrange(PATIENT_ID, SAMPLE_ID)

# Further limit the results by selecting adjacent ploidy states with the lowest clonality 
# values where multiple similar solutions are present. Plody values within a 
# threshold of 0.3 grouped together and the lowest clonality value is selected
pruned_results <- filtered_results %>%
  dplyr::arrange(SAMPLE_ID, ploidy) %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::mutate(pl_diff = abs(ploidy - dplyr::lag(ploidy))) %>% #, pu_diff = abs(purity - dplyr::lag(purity)) not used
  dplyr::mutate(new_state_n = dplyr::row_number() == 1 | pl_diff > 0.3) %>%
  dplyr::mutate(new_state = cumsum(new_state_n)) %>%
  dplyr::group_by(SAMPLE_ID, new_state) %>%
  dplyr::filter(rank_clonality == min(rank_clonality)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(use = rep(NA,times=nrow(.)),notes = rep(NA,times=nrow(.)))

write.table(filtered_results,paste0(outpath,project,"_filtered_results.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

write.table(pruned_results,paste0(outpath,project,"_fit_QC_predownsample.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

## ADDED by PS - adding output folder for results
if(!dir.exists(paste0(outpath,"plots"))){
	dir.create(paste0(outpath,"plots"))
}

if(length(unique(pruned_results$SAMPLE_ID)) == 1){
  i <- unique(pruned_results$SAMPLE_ID)
  dat <-  pruned_results %>%
    filter(SAMPLE_ID == i) %>%
    arrange(ploidy)
    #arrange(rank_clonality)
  x <- relative_smoothed
  cn <- assayDataElement(x,"copynumber")
  seg <- assayDataElement(x,"segmented")
  rel_ploidy <- mean(cn,na.rm=T)
  ll <- nrow(dat)
  
  png(paste0(outpath,"plots/", i, ".png"),type="cairo",w= 450*ll, h = 350)
  par(mfrow = c(1,ll))
  for(n in 1:nrow(dat)){
    ploidy <- dat[n,]$ploidy
    purity <- dat[n,]$purity
    cellploidy <- ploidy*purity+2*(1-purity)
    seqdepth <- rel_ploidy/cellploidy
    expTP53 <- round(dat[n,]$expected_TP53_AF, 2)
    TP53 <- dat[n,]$TP53freq

    #convert to abs

    pData(x)$ploidy <- ploidy
    pData(x)$purity <- purity

    temp <- x
    abs_cn <- depthtocn(cn,purity,seqdepth)
    abs_seg <- depthtocn(seg,purity,seqdepth)
    assayDataElement(temp,"copynumber") <- abs_cn
    assayDataElement(temp,"segmented") <- abs_seg

    #tmp_abs <- convert_rd_to_cn(x)
    # plot
    if(ploidy>5){
      yrange=15
    } else {
      yrange=10
    }
    plot(temp,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
           main=paste(i, " eTP53=",round(expTP53,2),
                      " AF=", round(TP53,2),
                      " p=",round(purity,2),
                      " pl=",round(ploidy,2),
                      sep=""),cex.main=0.8)
    abline(h=1:9,col = "blue")

  }
  dev.off()
} else {
  #relative_smoothed
  #Plot absolute CN fits for assessment -- currently multithreaded as all samples are processed together
  # Future implementation may change this to single thread and use scatter gather instead.
  foreach::foreach(i=unique(pruned_results$SAMPLE_ID)) %dopar%{
    dat <-  pruned_results %>%
      dplyr::filter(SAMPLE_ID == i) %>%
      dplyr::arrange(ploidy)
    x <- relative_smoothed[, i]
    cn <- Biobase::assayDataElement(x,"copynumber")
    seg <- Biobase::assayDataElement(x,"segmented")
    rel_ploidy <- mean(cn,na.rm=T)
    ll <- nrow(dat)
    png(paste0(outpath,"plots/", i, ".png"),type = "cairo", w= 450*ll, h = 350)
    par(mfrow = c(1,ll)) 
    for(n in 1:nrow(dat)){
      
      ploidy <- dat[n,]$ploidy
      purity <- dat[n,]$purity
      cellploidy <- ploidy * purity + (2*(1-purity))
      seqdepth <- rel_ploidy/cellploidy
      
      abs_cn <- depthtocn(cn,purity,seqdepth)
      abs_seg <- depthtocn(seg,purity,seqdepth)
      
      Biobase::assayDataElement(x,"copynumber") <- abs_cn
      Biobase::assayDataElement(x,"segmented") <- abs_seg
      
      expTP53 <- dat[n,]$expected_TP53_AF
      TP53 <- dat[n,]$TP53freq
      # plot   
      if(ploidy>5){
        yrange=15
      } else {
        yrange=10
      }
      plot(x,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
           main=paste(i, " eTP53=",round(expTP53,2),
                      " AF=", round(TP53,2),
                      " p=",round(purity,2),
                      " pl=",round(ploidy,2),
                      sep=""),cex.main=0.8)
      abline(h=1:yrange-1,col = "blue")
    }
    dev.off()
  }
}

#END
