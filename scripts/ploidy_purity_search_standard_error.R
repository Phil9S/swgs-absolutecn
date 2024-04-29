#this script outputs a ranked list of absolute copy number fits using the 
#traditional clonality error function. it is designed to be used in conjuction 
#with TP53 allele frequency to determine precise purity and ploidy fit

#grab commandline arguments
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[[1]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
cores <- as.numeric(snakemake@threads) 

metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")

# gridsearch params
pl_min <- snakemake@params[["ploidy_min"]] # default 1.6
pl_max <- snakemake@params[["ploidy_max"]] # default 8
pu_min <- snakemake@params[["purity_min"]] # default 0.15
pu_max <- snakemake@params[["purity_max"]] # deafult 1
#print(c(pl_min,pl_max,pu_min,pu_max))

#load libraries
suppressPackageStartupMessages(library(QDNAseqmod))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(matrixStats))
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))
registerDoMC(cores)
options(scipen = 999)

#read in relative copy number
rds.obj <- readRDS(rds.filename)
sampleNames(rds.obj)<-str_replace(sampleNames(rds.obj),"[(]","_")
sampleNames(rds.obj)<-str_replace(sampleNames(rds.obj),"[)]","_")
rds.pdata <- pData(rds.obj[[1]])

#set parameters for fixed bin size
bin_size <- bin*1000
bins<-getBinAnnotations(binSize = bin)
nbins_ref_genome <- sum(fData(rds.obj[[1]])$use)
nbins<-nrow(bins)

#define helper functions
depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

## TP53 target bin - hg19@30kb
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

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.obj[[1]])

#adjust for precomputed pl/pu
metaSample <- snakemake@wildcards[["sample"]]
metaflt <- metadata[metadata$SAMPLE_ID == metaSample,]

# set gridsearch limits based on availability of 
# precomputed pl/pu values
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
ploidies<-seq(pl_min,pl_max,0.1)
purities<-seq(pu_min,pu_max,0.01)

clonality<-c()
#ind<-which(colnames(rds.obj)==sample)
relcn<-rds.obj[[1]]
# added by PS
to_use <- fData(relcn)$use #
relcn <- relcn[to_use,] #
# added by PS
copynumber<-assayDataElement(relcn,"copynumber")
rel_ploidy<-mean(copynumber,na.rm=T)
num_reads<-sum(copynumber,na.rm=T)
sample <- sampleNames(relcn)

#print(sample)
res<-foreach(i=1:length(ploidies),.combine=rbind) %do% {
        ploidy<-ploidies[i]
	rowres<-foreach(j=1:length(purities),.combine=rbind) %do% {
            purity<-purities[j]
            downsample_depth<-(((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)
            cellploidy<-ploidy*purity+2*(1-purity)
            seqdepth<-rel_ploidy/cellploidy
            cn<-assayDataElement(relcn,"copynumber")
            seg<-assayDataElement(relcn,"segmented")
            cn<-as.numeric(cn[!is.na(cn),])
            seg<-as.numeric(seg[!is.na(seg),])
            integer_cn<-round(depthtocn(seg,purity,seqdepth))
            abs_cn<-depthtocn(seg,purity,seqdepth)
            diffs<-abs(abs_cn-integer_cn)
            TP53cn<-round(depthtocn(median(seg[gene_bin_seg]),purity,seqdepth),1) # to 1 decimal place
            expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))
            clonality<-mean(diffs)
            r <- c(ploidy,purity,clonality,downsample_depth,downsample_depth < rds.pdata$total.reads[row.names(rds.pdata)==sample],TP53cn,expected_TP53_AF)
	    r <- as.data.frame(t(r))
	    return(r)
        }
	return(as.data.frame(rowres))
}

colnames(res)<-c("ploidy","purity","clonality","downsample_depth","powered","TP53cn","expected_TP53_AF")
if(nrow(res) > 1){
	rownames(res)<-1:nrow(res)
	res<-data.frame(res,stringsAsFactors = F)
	res<-data.frame(apply(res,2,as.numeric,stringsAsFactors=F))
	res<-res[order(res$clonality,decreasing=FALSE),]
} else {
	rownames(res) <- 1
}

#output plot of clonality error landscape
pdf(snakemake@output[["pdf"]])
print(ggplot(res,aes(x=ploidy,y=purity,fill=clonality))+geom_tile()+
               scale_fill_gradient(low = "blue", high = "white",trans="log10")+
               theme_bw())
dev.off()

write.table(res,snakemake@output[["tsv"]],sep="\t",quote=F,row.names=FALSE)
