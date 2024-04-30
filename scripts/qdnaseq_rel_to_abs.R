# Clean env
args = commandArgs(trailingOnly=TRUE)

#load libraries
library(QDNAseqmod)
library(Biobase)
library(ggplot2)
library(stringr)
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))

qc.data <- read.table(snakemake@input[["meta"]],header = T,sep = "\t")
output_dir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
cores <- as.numeric(snakemake@threads)
registerDoMC(cores)

qc.data <- qc.data[qc.data$use == "TRUE",]
rds.filename <- snakemake@input[["rds"]]

rds.list <- lapply(rds.filename,FUN=function(x){readRDS(x)})

collapse_rds <- function(rds.list){
  comb <- rds.list[[1]][[1]]
  if(length(rds.list) > 1){
    for(i in 2:length(rds.list)){
      add <- rds.list[[i]][[1]]
      comb <- combine(comb,add)
    }
    rds.obj <- comb
  } else {
    rds.obj <- comb
  }
  return(rds.obj)
}

## TP53 target bin
target <- c("17:7565097-7590863")
get_gene_seg <- function(target=NULL,abs_data=NULL){
  to_use <- fData(abs_data)$use
  cn_obj <- abs_data[to_use,]
  bin_pos <- fData(cn_obj)[,c("chromosome","start","end")]
  chr <- as.numeric(str_split(string = target,pattern = ":",simplify = T)[1])
  start <- as.numeric(str_split(string = target,pattern = ":|-",simplify = T)[2])
  end <- as.numeric(str_split(string = target,pattern = ":|-",simplify = T)[3])

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
pData(rds.rel)$purity <- samples$purity[match(pData(rds.rel)$name,samples$SAMPLE_ID)]
pData(rds.rel)$ploidy <- samples$ploidy[match(pData(rds.rel)$name,samples$SAMPLE_ID)]
pData(rds.rel)$TP53freq <- samples$TP53freq[match(pData(rds.rel)$name,samples$SAMPLE_ID)]
pData(rds.rel)$PATIENT_ID <- samples$PATIENT_ID[match(pData(rds.rel)$name,samples$SAMPLE_ID)]

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.rel)

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
  TP53cn<-round(depthtocn(median(seg[gene_bin_seg]),purity,seqdepth),1) # to 1 decimal place / altered to correct bin value
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

# save segTable
getSegTable<-function(x){
    if(inherits(x,what = "QDNAseqCopyNumbers",which = F)){
        sn<-Biobase::assayDataElement(x,"segmented")
        fd <- Biobase::fData(x)
        fd$use -> use
        fdfiltfull<-fd[use,]
        sn<-sn[use,]
	      if(is.null(ncol(sn))){
		      sampleNa <- Biobase::sampleNames(x)
		      sn <- as.data.frame(sn)
		      colnames(sn) <- sampleNa
	      }
        segTable<-c()
        for(s in colnames(sn)){
            for(c in unique(fdfiltfull$chromosome))
            {
                snfilt<-sn[fdfiltfull$chromosome==c,colnames(sn) == s]
                fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
                sn.rle<-rle(snfilt)
                starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
                ends <- cumsum(sn.rle$lengths)
                lapply(1:length(sn.rle$lengths), function(s) {
                    from <- fdfilt$start[starts[s]]
                    to <- fdfilt$end[ends[s]]
                    segValue <- sn.rle$value[s]
                    c(fdfilt$chromosome[starts[s]], from, to, segValue)
                }) -> segtmp
                segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),sample = rep(s,times=nrow(matrix(unlist(segtmp), ncol=4, byrow=T))),stringsAsFactors=F)
                segTable<-rbind(segTable,segTableRaw)
            }
        }
        colnames(segTable) <- c("chromosome", "start", "end", "segVal","sample")
        return(segTable)
    } else {
        # NON QDNASEQ BINNED DATA FUNCTION
	stop("segtable error")
    }
}

write.table(getSegTable(abs_profiles),
	paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/",project,"_",bin,"kb_ds_absCopyNumber_segTable.tsv"),
	sep = "\t",quote=F,row.names=FALSE)

#write table of fits
write.table(res,paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/",project,"_",bin,"kb_ds_abs_fits.tsv"),sep = "\t",quote=F,row.names=FALSE)
