# gridsearch.R 
## Outputs a table of absolute copy number fits across a ploidy and purity gridsearch
## using mean absolute error (aka clonality) error function, sorted by lowest.
## it is designed to be used in conjunction with TP53 allele frequency to determine
## precise purity and ploidy fit.

# grab commandline arguments passed via snakemake object
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[["rds"]]
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
genome <- as.character(snakemake@params[["genome"]]) # hg19 or hg38
pl_min <- snakemake@params[["ploidy_min"]] # default 1.6 
pl_max <- snakemake@params[["ploidy_max"]] # default 8
pu_min <- snakemake@params[["purity_min"]] # default 0.15
pu_max <- snakemake@params[["purity_max"]] # default 1
sampleName <- as.character(snakemake@params[["sample"]])
af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
filter_underpowered <- as.logical(snakemake@params[["filter_underpowered"]]) # Filter for powered fits only
filter_homozygous <- as.logical(snakemake@params[["filter_homozygous"]]) # filter homozygous loss
homozygous_prop <- as.numeric(snakemake@params[["homozygous_prop"]])
homozygous_thrsh <- snakemake@params[["homozygous_threshold"]] # default 0.4

outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_PRE_down_sampling/")

# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)
# read in relative copy number and extract model and read data
rds.obj <- readRDS(rds.filename) 
total_reads <- Biobase::pData(rds.obj)$total.reads
# set whitelisted ref_genome_bins for fixed bin size
bin_size <- bin * 1000
nbins_ref_genome <- sum(Biobase::fData(rds.obj)$use)
#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(abs_data = rds.obj,genome = genome)
# set gridsearch limits based on availability of precomputed pl/pu values
metaflt <- metadata[metadata$SAMPLE_ID == sampleName,]

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
ploidies <- seq.int(pl_min,pl_max,0.1)
purities <- seq.int(pu_min,pu_max,0.01)

# Extract CN and Seg data
to_use <- Biobase::fData(rds.obj)$use
relcn <- rds.obj[to_use,]

# Iterate through ploidy/purity combinations
res <- data.frame()
for(i in 1:length(ploidies)){
  ploidy <- ploidies[i]
  rowres <- data.frame()
	for(j in 1:length(purities)){
    purity <- purities[j]
    
    gridVals <- gridStats(obj = relcn,ploidy = ploidy,purity = purity)
    
    downsample_depth <- getDownsampleDepth(ploidy,purity,nbins_ref_genome)
    num_segs <- length(rle(as.numeric(gridVals$seg))$values)
    
    errors <- getErrors(data = gridVals)
    
    targetCNVal <- median(gridVals$seg[gene_bin_seg])
    TP53cn <- round(depthtocn(targetCNVal,purity,gridVals$seqdepth),1) # to 1 decimal place
    expected_TP53_AF <- TP53cn * purity / (TP53cn * purity + 2 * (1-purity))
    
    hmzyg <- sum(gridVals$abs_seg <= homozygous_thrsh) * bin_size
    powered <- downsample_depth < total_reads
    
    r <- c(ploidy,purity,num_segs,errors,downsample_depth,
           powered,TP53cn,expected_TP53_AF,hmzyg)
    r <- as.data.frame(t(r))
    rowres <- rbind(rowres,r)
	}
  res <- rbind(res,rowres)
}
# Format gridsearch table
res <- cbind(rep(metaSample,times=nrow(res)),res)

colnames(res) <- fittingColumnNames
res <- res[order(res$clonality,decreasing=FALSE),]

# Output sunrise plot of clonality error landscape
pdf(snakemake@output[["pdf"]])
print(ggplot2::ggplot(res,ggplot2::aes(x=ploidy,y=purity,fill=clonality)) +
        ggplot2::geom_tile()+
        ggplot2::scale_x_continuous(expand = c(0,0),breaks = ploidies[seq.int(1,length(ploidies),4)]) +
        ggplot2::scale_y_continuous(expand = c(0,0),breaks = purities[seq.int(1,length(purities),5)]) +
        ggplot2::scale_fill_gradient(low = "blue", high = "white",name = "clonality\n(MAE)") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom"))
dev.off()

write.table(res,snakemake@output[["tsv"]],sep="\t",quote=F,row.names=FALSE)

### Gridsearch filtering
filteredTables <- filterFitTable(table = res,metadata = metadata,
                                 filter_underpowered = filter_underpowered,
                                 filter_homozygous = filter_homozygous,
                                 af_cutoff = af_cutoff)

write.table(filteredTables$filtered,snakemake@output[["filt"]],
            sep="\t",col.names=T,row.names=F,quote=F)

write.table(filteredTables$pruned,snakemake@output[["fit"]],
            sep="\t",col.names=T,row.names=F,quote=F)

## Plot fits
if(!dir.exists(paste0(outpath,"plots"))){
  dir.create(paste0(outpath,"plots"),recursive = TRUE)
}

ll <- nrow(filteredTables$pruned)
png(paste0(outpath,"plots/",sampleName,".png"),type = "cairo", w= 450*ll, h = 350)
par(mfrow = c(1,ll)) 
for(n in 1:nrow(pruned_results)){
  
  ploidy <- filteredTables$pruned[n,]$ploidy
  purity <- filteredTables$pruned[n,]$purity
  
  gridValsPlot <- gridStats(obj = relcn,ploidy = ploidy,purity = purity)
  errors <- getErrors(gridValsPlot)
  
  Biobase::assayDataElement(relcn,"copynumber") <- gridValsPlot$abs_cn
  Biobase::assayDataElement(relcn,"segmented") <- gridValsPlot$abs_seg
  
  plotProfile(relcn,ploidy = ploidy,purity = purity,
              clonality = errors["clonality"],rmse = errors["rmse"])
}
dev.off()

#END
