# QDNAseqmod relative to absolute fitting
args = commandArgs(trailingOnly=TRUE)

rds.filename <- snakemake@input[["rds"]]
metafile <- snakemake@input[["meta"]]
metadata <- read.table(metafile,header = T,sep = "\t")
output_dir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
genome <- snakemake@params[["genome"]]

source("scripts/funcs.R")

outpath <- paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/")
metadata <- metadata[metadata$use == "TRUE",]

rds.list <- lapply(rds.filename,FUN = function(x){
  readRDS(x)
})

## TP53 target bin (genome dependent)
if(genome == "hg19"){
	target <- c("17:7565097-7590863")
} else if(genome == "hg38"){
	target <- c("17:7661779-7687538")
}

# Combine and load rds objects
rds.rel <- collapse_rds(rds.list)

if(!dir.exists(paste0(outpath,"plots"))){
	dir.create(paste0(outpath,"plots"))
}

# List samples
samples <- metadata[which(metadata$SAMPLE_ID %in% colnames(rds.rel)),]

# Add pheno information
Biobase::pData(rds.rel)$purity <- samples$purity[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]
Biobase::pData(rds.rel)$ploidy <- samples$ploidy[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]
Biobase::pData(rds.rel)$TP53freq <- samples$TP53freq[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]
Biobase::pData(rds.rel)$PATIENT_ID <- samples$PATIENT_ID[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]
Biobase::pData(rds.rel)$rmse <- samples$rmse[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]
Biobase::pData(rds.rel)$mae <- samples$clonality[match(Biobase::pData(rds.rel)$name,samples$SAMPLE_ID)]

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.rel)

# Generate abs plot and table of fits
res <- data.frame()
abs_profiles <- rds.rel[Biobase::fData(rds.rel)$use,]
# For each
for(sample in Biobase::pData(abs_profiles)$name){
  # Index and subselect sample
  #ind <- which(colnames(rds.rel) == sample)
  relcn <- abs_profiles[,sample]
  #to_use <- fData(relcn)$use
  #relcn <- relcn[to_use,]
  #smooth.bool <- FALSE
  # Extract cn and ploidy
  
  # Extract CN and Segs
  cn <- Biobase::assayDataElement(relcn,"copynumber")
  seg <- Biobase::assayDataElement(relcn,"segmented")
  
  rel_ploidy <- mean(cn,na.rm=T)
  ploidy <- Biobase::pData(relcn)$ploidy
  purity <- Biobase::pData(relcn)$purity
  cellploidy <- ploidy * purity + 2*(1 - purity)
  seqdepth <- rel_ploidy / cellploidy
  
  
  # Convert to abs
  abs_cn <- depthtocn(cn,purity,seqdepth)
  abs_seg <- depthtocn(seg,purity,seqdepth)
  
  Biobase::assayDataElement(relcn,"copynumber") <- abs_cn
  Biobase::assayDataElement(relcn,"segmented") <- abs_seg
  
  # Add TP53 info
  targetCNVal <- median(seg[gene_bin_seg])
  TP53cn <- round(depthtocn(targetCNVal,purity,seqdepth),1) # to 1 decimal place / altered to correct bin value
  expected_TP53_AF <- TP53cn * purity / (TP53cn * purity + 2*(1 - purity))
  TP53freq <- Biobase::pData(relcn)$TP53freq
  
  # Add patient-level info
  pat <- as.character(Biobase::pData(relcn)$PATIENT_ID)
  res <- rbind(res,matrix(c(sample,pat,ploidy,purity,TP53cn,
                            round(expected_TP53_AF,2),TP53freq,NA,NA),nrow = 1,ncol = 9))
  
  # Y axis range
  if(ploidy > 5){
    yrange = 15
  } else {
    yrange = 10
  }
  # Plot abs fit
  
  mae <- Biobase::pData(relcn)$mae
  rmse <- Biobase::pData(relcn)$rmse
  sub <- paste0("ploidy=",round(ploidy,2)," | ",
                " purity=",round(purity,2)," | ",
                " MAE=",round(mae,3)," | ",
                " RMSE=",round(rmse,3))
  
  png(paste0(outpath,"plots/",sample,".png"),type="cairo",w = 8,h = 6,unit="in",res = 250)
  par(mfrow = c(1,1))
  plot(relcn,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,
       ylim=c(0,yrange),ylab="Absolute tumour CN")
  abline(h=1:yrange-1, col = "blue")
  mtext(sub,side = 1,line = 3.5)
  dev.off()
  
  # Add to abs RDS
  Biobase::assayDataElement(abs_profiles,"copynumber")[,sample] <- abs_cn
  Biobase::assayDataElement(abs_profiles,"segmented")[,sample] <- abs_seg
}

# Annotated and rename table
colnames(res) <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","TP53cn","expected_TP53_AF","TP53freq","use","notes")
res <- data.frame(res,stringsAsFactors = F)

# Save rds
saveRDS(abs_profiles,file=paste0(outpath,project,"_",bin,"kb_ds_absCopyNumber.rds"))

# save segTable
segTable <- getSegTable(abs_profiles)
write.table(segTable,paste0(outpath,project,"_",bin,"kb_ds_absCopyNumber_segTable.tsv"),
	sep = "\t",quote=F,row.names=FALSE)

#write table of fits
write.table(res,paste0(outpath,,project,"_",bin,"kb_ds_abs_fits.tsv"),
  sep = "\t",quote=F,row.names=FALSE)
