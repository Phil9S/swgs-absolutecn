library(tidyverse)
library(Biobase)
library(QDNAseqmod)
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))

## Added by PS
args = commandArgs(trailingOnly=TRUE)

metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]

#threads
cores <- as.numeric(snakemake@threads)
registerDoMC(cores)

# Filter parameters
filter_underpowered <- snakemake@params[["filter_underpowered"]]
af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
filter_homozygous=snakemake@params[["filter_homozygous"]]
homozygous_prop=as.numeric(snakemake@params[["homozygous_prop"]])

# filter homozygous loss
if(filter_homozygous == "TRUE"){
  filter_homozygous <- TRUE
} else {
  filter_homozygous <- FALSE
} 

# Filter for powered fits only
if(filter_underpowered == "TRUE"){
  filter_underpowered <- TRUE
} else {
  filter_underpowered <- FALSE
}

# read in relative CN data
# collapse rds files function
rds.filename <- snakemake@input[["rds"]]

rds.list <- lapply(rds.filename,FUN=function(x){readRDS(x)})
collapse_rds <- function(rds.list){
  comb <- rds.list[[1]][[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]][[1]]
      comb <- Biobase::combine(comb,add)
    }
    rds.obj <- comb
  } else {
    rds.obj <- comb
  } 
  return(rds.obj)
}
# Combine and load rds objects
relative_smoothed <- collapse_rds(rds.list)
saveRDS(relative_smoothed,paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_",bin,"kb_relSmoothedCN.rds"))

filelist <- snakemake@input[["cl"]]
clonality <- do.call(rbind,
			lapply(filelist,FUN = function(x){
				n <- gsub(pattern="_clonality.tsv",rep="",x=x)
        prefix <- paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",project,"_")
        n <- gsub(pattern=prefix,rep="",x=n)
				tab <- read.table(x,sep="\t",skip=1)
				tab <- cbind(rep(n,times=nrow(tab)),tab)
				return(tab)
			}))
colnames(clonality) <- c("SAMPLE_ID","ploidy","purity","clonality","downsample_depth","powered","TP53cn","expected_TP53_AF","homozygousLoss")

clonality <- left_join(clonality,metadata,by="SAMPLE_ID") %>%
                select(SAMPLE_ID,PATIENT_ID,ploidy,purity,clonality,downsample_depth,powered,TP53cn,expected_TP53_AF,TP53freq,smooth,homozygousLoss)

depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
  (x/seqdepth-2*(1-purity))/purity
}

#Top 10 when ranking by clonality and TP53
filtered_results <- clonality

## Filter underpowered
if(filter_underpowered){
  filtered_results <- filtered_results %>%
    filter(powered == 1) 
}


## filter homozygous loss
if(filter_homozygous){
  filtered_results <- filtered_results %>%
    filter(homozygousLoss <= homozygous_prop)
}

# standard filtering
filtered_results <- filtered_results %>% # filter underpowered fits when config variable is TRUE
  select(SAMPLE_ID, PATIENT_ID, everything()) %>%
  group_by(SAMPLE_ID, ploidy) %>%
  mutate(rank_clonality = min_rank(clonality)) %>% #rank clonality within a unique ploidy state 
  filter(rank_clonality == 1) %>% #select ploidy with the lowest clonality within a unique ploidy state 
  group_by(SAMPLE_ID) %>%
  top_n(-10, wt = clonality) %>% # select top 10 ploidy states with the lowest clonality values
  mutate(rank_clonality = min_rank(clonality)) %>% # rank by clonality within a sample across ploidies in top 10
  # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
  filter(is.na(TP53freq) | near(expected_TP53_AF,TP53freq, tol = af_cutoff )) %>%  
  arrange(PATIENT_ID, SAMPLE_ID)

#Further limit the results by selecting the ploidy states with the lowest clonality values where multiple similar solutions are present.
#Threshold of 0.3 used to select different states

pruned_results <- filtered_results %>%
  arrange(SAMPLE_ID, ploidy) %>%
  group_by(SAMPLE_ID) %>%
  mutate(pl_diff = abs(ploidy-lag(ploidy)), pu_diff = abs(purity-lag(purity))) %>%
  mutate(new = row_number() == 1 | pl_diff > 0.3) %>%
  mutate(new_state = cumsum(new)) %>%
  group_by(SAMPLE_ID, new_state) %>%
  filter(rank_clonality == min(rank_clonality))

pruned_results$use <- rep(NA,times=nrow(pruned_results))
pruned_results$notes <- rep(NA,times=nrow(pruned_results))

write.table(filtered_results,
  paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_filtered_results.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

write.table(pruned_results,
  paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/",project,"_fit_QC_predownsample.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

## ADDED by PS - adding output folder for results
if(!dir.exists(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots"))){
	dir.create(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots"))
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
  png(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots/", i, ".png"),
	type="cairo",w= 450*ll, h = 350)
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
  #Plot absolute CN fits for assessment
  foreach(i=unique(pruned_results$SAMPLE_ID)) %dopar% {
    dat <-  pruned_results %>%
      filter(SAMPLE_ID == i) %>%
      arrange(ploidy)
      #arrange(rank_clonality)
    x <- relative_smoothed[, i]
    cn <- assayDataElement(x,"copynumber")
    seg <- assayDataElement(x,"segmented")
    rel_ploidy <- mean(cn,na.rm=T)
    ll <- nrow(dat)
    png(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/plots/", i, ".png"),
	type = "cairo", w= 450*ll, h = 350,)
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
    abline(h=1:yrange-1,col = "blue")
  
    }
    dev.off()
  }
}

#END
