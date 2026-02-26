# gridsearch.R 
## Outputs a table of absolute copy number fits across a ploidy and purity gridsearch
## using mean absolute error (aka clonality) error function, sorted by lowest.
## it is designed to be used in conjunction with TP53 allele frequency to determine
## precise purity and ploidy fit

# grab commandline arguments passed via snakemake object
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[[1]]
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
genome <- as.character(snakemake@params[["genome"]]) # hg19 or hg38
hmz_thrsh <- snakemake@params[["homozygous_threshold"]] # default 0.4
pl_min <- snakemake@params[["ploidy_min"]] # default 1.6 
pl_max <- snakemake@params[["ploidy_max"]] # default 8
pu_min <- snakemake@params[["purity_min"]] # default 0.15
pu_max <- snakemake@params[["purity_max"]] # default 1

# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)

rds.obj <- readRDS(rds.filename) # read in relative copy number and extract model and read data
total_reads <- Biobase::pData(rds.obj)$total.reads

# set whitelisted ref_genome_bins for fixed bin size
bin_size <- bin * 1000
nbins_ref_genome <- sum(Biobase::fData(rds.obj)$use)

## TP53 target bin (genome dependent) hard coded - add as parameter?
if(genome == "hg19"){
        target <- c("17:7565097-7590863")
} else if(genome == "hg38"){
        target <- c("17:7661779-7687538")
}
#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.obj)

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
    
    r <- c(ploidy,purity,clonality,rmse,downsample_depth,powered,TP53cn,expected_TP53_AF,hmzyg,MedianSegVar)
    r <- as.data.frame(t(r))
    rowres <- rbind(rowres,r)
	}
  res <- rbind(res,rowres)
}
# Format gridsearch table
colnames(res) <- c("ploidy","purity","clonality","rmse","downsample_depth",
                   "powered","TP53cn","expected_TP53_AF","homozygousLoss","MedianSegVar")
res <- res[order(res$clonality,decreasing=FALSE),]

# Output sunrise plot of clonality error landscape
pdf(snakemake@output[["pdf"]])
print(ggplot2::ggplot(res,ggplot2::aes(x=ploidy,y=purity,fill=clonality)) +
        ggplot2::geom_tile()+
        ggplot2::scale_x_continuous(expand = c(0,0),breaks = ploidies[seq.int(1,length(ploidies),2)]) +
        ggplot2::scale_y_continuous(expand = c(0,0),breaks = purities[seq.int(1,length(purities),5)]) +
        ggplot2::scale_fill_gradient(low = "blue", high = "white",name = "clonality\n(MAE)") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom"))
dev.off()

write.table(res,snakemake@output[["tsv"]],sep="\t",quote=F,row.names=FALSE)
# END