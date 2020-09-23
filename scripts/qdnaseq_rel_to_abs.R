# Clean env
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

#load libraries
library(QDNAseqmod)
library(Biobase)
library(ggplot2)
library(stringr)
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))

qc.data <- read.table(args[1],header = T,sep = "\t")
output_dir <- args[2]
bin <- as.numeric(args[3])
project <- args[4]
cores <- args[5]

registerDoMC(cores)

qc.data <- qc.data[qc.data$use == "TRUE",]

rds.filename <- list.files(pattern="*relSmoothedCN.rds",path=paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/relative_cn_rds/"))
rds.list <- lapply(rds.filename,FUN=function(x){readRDS(paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/relative_cn_rds/",x))})

collapse_rds <- function(rds.list){
  comb <- rds.list[[1]][[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]][[1]]
      comb <- combine(comb,add)
    }
    rds.obj <- comb
  }
  return(rds.obj)
}

# Combine and load rds objects
rds.rel <- collapse_rds(rds.list)

if(!dir.exists(paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/plots"))){
	dir.create(paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/plots"))
}

# convert depth to abs cn
depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
  (x/seqdepth-2*(1-purity))/purity
}

# List samples
samples <- qc.data[which(qc.data$SAMPLE_ID %in% colnames(rds.rel)),]

# Add pheno information
pData(rds.rel)$purity <- samples$purity[match(samples$SAMPLE_ID,pData(rds.rel)$name)]
pData(rds.rel)$ploidy <- samples$ploidy[match(samples$SAMPLE_ID,pData(rds.rel)$name)]
pData(rds.rel)$TP53freq <- samples$TP53freq[match(samples$SAMPLE_ID,pData(rds.rel)$name)]
pData(rds.rel)$PATIENT_ID <- samples$PATIENT_ID[match(samples$SAMPLE_ID,pData(rds.rel)$name)]

# Generate abs plot and table of fits
res <- data.frame(matrix(ncol = 9, nrow = 0))
abs_profiles <- rds.rel[fData(rds.rel)$use,]
# For each
for(sample in pData(rds.rel)$name){
  # Index and subselect sample
  ind <- which(colnames(rds.rel)==sample)
  relcn <- rds.rel[,ind]
  to_use <- fData(relcn)$use #
  relcn <- relcn[to_use,]
  smooth.bool <- FALSE
  # Extract cn and ploidy
  copynumber <- assayDataElement(relcn,"copynumber")
  rel_ploidy <- mean(copynumber,na.rm=T)
  ploidy <- pData(relcn)$ploidy
  purity <- pData(relcn)$purity
  cellploidy <- ploidy*purity+2*(1-purity)
  seqdepth <- rel_ploidy/cellploidy

  # Extract CN and Segs
  cn <- assayDataElement(relcn,"copynumber")
  seg <- assayDataElement(relcn,"segmented")
  
  # Convert to abs
  abs_cn <- depthtocn(cn,purity,seqdepth)
  abs_seg <- depthtocn(seg,purity,seqdepth)
  assayDataElement(relcn,"copynumber") <- abs_cn
  assayDataElement(relcn,"segmented") <- abs_seg
  # Add to abs RDS
  assayDataElement(abs_profiles,"copynumber")[,ind] <- abs_cn
  assayDataElement(abs_profiles,"segmented")[,ind] <- abs_seg
  # Add TP53 info
  TP53cn<-round(depthtocn(seg[73504],purity,seqdepth),1) # to 1 decimal place / altered to correct bin value
  expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))
  TP53freq <- pData(relcn)$TP53freq
  # Add patient-level info
  pat <- as.character(pData(relcn)$PATIENT_ID)
  res <- rbind(res,matrix(c(sample,pat,ploidy,purity,TP53cn,round(expected_TP53_AF,2),TP53freq,NA,NA),nrow = 1,ncol = 9))
  
  # Y axis range
  if(ploidy>5){
    yrange=15
  } else {
    yrange=10
  }
  # Plot abs fit
  png(paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/plots/",sample,".png"), w= 8, h = 6, unit="in", res = 250)
  par(mfrow = c(1,1))
  plot(relcn,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
       main=paste(sample, " eTP53=",round(expected_TP53_AF,2),
                  " AF=", round(TP53freq,2),
                  " p=",round(purity,2),
                  " pl=",round(ploidy,2),
                  sep=""))
  abline(h=1:9, col = "blue")
  dev.off()
}

# Annotated and rename table
colnames(res) <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","TP53cn","expected_TP53_AF","TP53freq","use","notes")
res <- data.frame(res,stringsAsFactors = F)

# Save rds
saveRDS(abs_profiles,file=paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/",project,"_",bin,"kb_ds_absCopyNumber.rds"))

#write table of fits
write.table(res,paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/",project,"_",bin,"kb_ds_abs_fits.tsv"),sep = "\t",quote=F,row.names=FALSE)
