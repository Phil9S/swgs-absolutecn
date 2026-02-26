#this script outputs a ranked list of absolute copy number fits using the 
#traditional clonality error function. it is designed to be used in conjuction 
#with TP53 allele frequency to determine precise purity and ploidy fit

#grab commandline arguments
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[[1]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
genome <- as.character(snakemake@params[["genome"]])
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")

# gridsearch params
pl_min <- snakemake@params[["ploidy_min"]] # default 1.6
pl_max <- snakemake@params[["ploidy_max"]] # default 8
pu_min <- snakemake@params[["purity_min"]] # default 0.15
pu_max <- snakemake@params[["purity_max"]] # default 1

# homozygous threshold
hmz_thrsh <- snakemake@params[["homozygous_threshold"]]

#load libraries
suppressPackageStartupMessages(library(QDNAseqmod))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(matrixStats))
options(scipen = 999)

### test file 
rds.filename <- "../../../Downloads/qdnaseqmod_update_pseudobulk_B4_30kb_relSmoothedCN.rds"
###

# read in relative copy number
rds.obj <- readRDS(rds.filename)
# renaming duplicated and potentailly not needed
# sampleNames(rds.obj) <- str_replace(sampleNames(rds.obj),"[(]","_")
rds.pdata <- Biobase::pData(rds.obj)
total_reads <- rds.pdata$total.reads

# set whitelisted ref_genome_bins for fixed bin size
bin_size <- bin * 1000
nbins_ref_genome <- sum(Biobase::fData(rds.obj)$use)

# unused
#bins <- getBinAnnotations(binSize = bin)
#nbins <- nrow(bins) unused

#define helper functions
# converts readdepth to copy number given purity and single copy depth
depthtocn<-function(x,purity,seqdepth){
    (x/seqdepth-2*(1-purity))/purity
}

get_gene_seg <- function(target=NULL,abs_data=NULL){
  to_use <- Biobase::fData(abs_data)$use
  cn_obj <- abs_data[to_use,]
  bin_pos <- Biobase::fData(cn_obj)[,c("chromosome","start","end")]
  #chr <- as.numeric(stringr::str_split(string = target,pattern = ":",simplify = T)[1])
  #start <- as.numeric(stringr::str_split(string = target,pattern = ":|-",simplify = T)[2])
  #end <- as.numeric(stringr::str_split(string = target,pattern = ":|-",simplify = T)[3])
  # Adjust target splitting and assign variables using mapply
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

## TP53 target bin (genome dependent)
## hard coded - add as parameter?
if(genome == "hg19"){
        target <- c("17:7565097-7590863")
} else if(genome == "hg38"){
        target <- c("17:7661779-7687538")
}

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.obj)

# adjust for precomputed pl/pu
# set gridsearch limits based on availability of precomputed pl/pu values
metaSample <- Biobase::sampleNames(rds.obj)
metaflt <- metadata[metadata$SAMPLE_ID == metaSample,]

if(!is.null(metaflt$precPloidy)){
	if(!is.na(metaflt$precPloidy)){
    pl_min <- metaflt$precPloidy
		pl_max <- metaflt$precPloidy
  }
}

if(!is.null(metaflt$precPurity)){
	if(!is.na(metaflt$precPurity)){
    pu_min <- metaflt$precPurity
		pu_max <- metaflt$precPurity
	}
}

#estimate absolute copy number fits for all samples in parallel
ploidies <- seq.int(pl_min,pl_max,0.1)
purities <- seq.int(pu_min,pu_max,0.01)

#clonality<-c()
#relcn <- rds.obj
# Extract CN and Seg data
to_use <- Biobase::fData(rds.obj)$use
relcn <- rds.obj[to_use,]
seg <- Biobase::assayDataElement(relcn,"segmented")
seg <- as.numeric(seg[!is.na(seg),])
cn <- Biobase::assayDataElement(relcn,"copynumber")
cn <- as.numeric(cn[!is.na(cn),])

# Compute readcount ploidy and read counts
rel_ploidy <- mean(cn,na.rm=T)
num_reads <- sum(cn,na.rm=T)
targetCNVal <- median(seg[gene_bin_seg])

res <- data.frame()
for(i in 1:length(ploidies)){
  ploidy <- ploidies[i]
  rowres <- data.frame()
	for(j in 1:length(purities)){
    purity<-purities[j]
    downsample_depth <- getDownsampleDepth(ploidy,purity,nbins_ref_genome)
    print(downsample_depth)
    cellploidy <- purity * ploidy + (2 * (1 - purity))
    seqdepth <- rel_ploidy / cellploidy
    
    abs_cn <- depthtocn(seg,purity,seqdepth)
    abs_cnbin <- depthtocn(cn,purity,seqdepth)
    integer_cn <- round(abs_cn,digits = 0)
    
    errors <- abs_cn - integer_cn
    
    TP53cn <- round(depthtocn(targetCNVal,purity,seqdepth),1) # to 1 decimal place
    expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))
    
    clonality <- mean(abs(errors)) # clonality is a legacy name for MAE
    rmse <- sqrt(mean(errors^2)) # Root Mean Squared Error
    MedianSegVar <- calculateSegmentVar(abs_cn = abs_cn,abs_cnbin = abs_cnbin)
    
    hmzyg <- sum(abs_cn <= hmz_thrsh) * bin_size
    powered <- downsample_depth < total_reads
    
    r <- c(ploidy,purity,clonality,rmse,downsample_depth,powered,TP53cn,expected_TP53_AF,hmzyg)
    r <- as.data.frame(t(r))
    rowres <- rbind(rowres,r)
	}
  res <- rbind(res,rowres)
}

colnames(res) <- c("ploidy","purity","clonality","rmse","downsample_depth",
                   "powered","TP53cn","expected_TP53_AF","homozygousLoss")

# These data.frame operations don't appear to do anything
#rownames(res) <- 1:nrow(res)
#res<-data.frame(res,stringsAsFactors = F)
#res<-data.frame(apply(res,2,as.numeric,stringsAsFactors=F))
res <- res[order(res$clonality,decreasing=FALSE),]

#output plot of clonality error landscape
pdf(snakemake@output[["pdf"]])
print(ggplot2::ggplot(res,ggplot2::aes(x=ploidy,y=purity,fill=clonality)) +
        ggplot2::geom_tile()+
        ggplot2::scale_x_continuous(expand = c(0,0),breaks = ploidies[seq.int(1,length(ploidies),2)]) +
        ggplot2::scale_y_continuous(expand = c(0,0),breaks = purities[seq.int(1,length(purities),5)]) +
        ggplot2::scale_fill_gradient(low = "blue", high = "white",name = "clonality\n(MAE)") +
        ggplot2::theme_bw())
dev.off()

write.table(res,snakemake@output[["tsv"]],sep="\t",quote=F,row.names=FALSE)
