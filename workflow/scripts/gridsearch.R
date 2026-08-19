# gridsearch.R 
## Outputs a table of absolute copy number fits across a ploidy and purity gridsearch
## using mean absolute error (aka clonality) error function, sorted by lowest.
## it is designed to be used in conjunction with TP53 allele frequency to determine
## precise purity and ploidy fit.

# grab commandline arguments passed via snakemake object
args = commandArgs(trailingOnly=TRUE)

if(exists("snakemake")){
  rds <- snakemake@input[["rds"]]
  meta <- snakemake@params[["meta"]]
  
  pdfOut <- snakemake@output[["pdf"]]
  tsvOut <- snakemake@output[["tsv"]]
  plotOut <- snakemake@output[["plot"]]
  filtOut <- snakemake@output[["filt"]]
  fitOut <- snakemake@output[["fit"]]
  
  bin <- as.numeric(snakemake@params[["bin"]])
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
  homozygous_thrsh <- as.numeric(snakemake@params[["homozygous_threshold"]]) # default 0.4
  
} else {
  
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "gridsearch"))
  
  rds <- opts$rds
  meta <- opts$meta
  
  pdfOut <- opts$pdfOut
  tsvOut <- opts$tsvOut
  plotOut <- opts$plotOut
  filtOut <- opts$filtOut
  fitOut <- opts$fitOut
  
  bin <- as.numeric(opts$bin)
  project <- opts$project
  genome <- opts$genome
  
  pl_min <- opts$pl_min
  pl_max <- opts$pl_max
  pu_min <-opts$pu_min
  pu_max <- opts$pu_max
  sampleName <- as.character(opts$sampleName)
  
  af_cutoff <- as.numeric(opts$af_cutoff)
  
  filter_underpowered <- as.logical(opts$filter_underpowered)
  filter_homozygous <- as.logical(opts$filter_homozygous)
  homozygous_prop <- as.numeric(opts$homozygous_prop)
  homozygous_thrsh <- as.numeric(opts$homozygous_thrsh)
}

metadata <- read.table(file = meta,header=T,sep="\t")

# source functions (temp until moved to package) and set scipen
options(scipen = 999)

# read in relative copy number and extract model and read data
rdsObj <- readRDS(rds) 
total_reads <- Biobase::pData(rdsObj)$total.reads

# set whitelisted ref_genome_bins for fixed bin size
bin_size <- bin * 1000
nbins_ref_genome <- sum(Biobase::fData(rdsObj)$use)

#Get target anchor gene segments
gene_bin_seg <- rswgsabsolutecn::getGeneSeg(abs_data = rdsObj,genome = genome)

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
to_use <- Biobase::fData(rdsObj)$use
relcn <- rdsObj[to_use,]

# Iterate through ploidy/purity combinations
res <- data.frame()
for(i in 1:length(ploidies)){
  ploidy <- ploidies[i]
  rowres <- data.frame()
  for(j in 1:length(purities)){
    purity <- purities[j]
    
    gridVals <- rswgsabsolutecn::gridStats(obj = relcn,ploidy = ploidy,purity = purity)
    
    downsample_depth <- rswgsabsolutecn::getDownsampleDepth(ploidy,purity,nbins_ref_genome)
    num_segs <- length(rle(as.numeric(gridVals$seg))$values)
    
    errors <- rswgsabsolutecn::getErrors(data = gridVals)
    
    targetCNVal <- median(gridVals$seg[gene_bin_seg])
    TP53cn <- round(rswgsabsolutecn::depthtocn(targetCNVal,purity,gridVals$seqdepth),1) # to 1 decimal place
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
res <- cbind(rep(sampleName,times=nrow(res)),res)

colnames(res) <- get(utils::data("fittingColumnNames",package = "rswgsabsolutecn",envir = environment()))
res <- res[order(res$clonality,decreasing=FALSE),]

# Output sunrise plot of clonality error landscape
pdf(pdfOut)
print(ggplot2::ggplot(res,ggplot2::aes(x=ploidy,y=purity,fill=clonality)) +
        ggplot2::geom_tile()+
        ggplot2::scale_x_continuous(expand = c(0,0),breaks = ploidies[seq.int(1,length(ploidies),4)]) +
        ggplot2::scale_y_continuous(expand = c(0,0),breaks = purities[seq.int(1,length(purities),5)]) +
        ggplot2::scale_fill_gradient(low = "blue", high = "white",name = "clonality\n(MAE)") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom"))
dev.off()

write.table(res,tsvOut,sep="\t",quote=F,row.names=FALSE)

### Gridsearch filtering
filteredTables <- rswgsabsolutecn::filterFitTable(table = res,metadata = metadata,
                                 filter_underpowered = filter_underpowered,
                                 filter_homozygous = filter_homozygous,
                                 homozygous_prop = homozygous_prop,af_cutoff = af_cutoff)

write.table(filteredTables$filtered,filtOut,
            sep="\t",col.names=T,row.names=F,quote=F)

write.table(filteredTables$pruned,fitOut,
            sep="\t",col.names=T,row.names=F,quote=F)

ll <- nrow(filteredTables$pruned)
png(plotOut,type = "cairo", w= 450*ll, h = 350)
par(mfrow = c(1,ll)) 
for(n in 1:nrow(filteredTables$pruned)){
  
  ploidy <- filteredTables$pruned[n,]$ploidy
  purity <- filteredTables$pruned[n,]$purity
  
  gridValsPlot <- rswgsabsolutecn::gridStats(obj = relcn,ploidy = ploidy,purity = purity)
  errors <- rswgsabsolutecn::getErrors(gridValsPlot)
  
  Biobase::assayDataElement(relcn,"copynumber") <- gridValsPlot$abs_cn
  Biobase::assayDataElement(relcn,"segmented") <- gridValsPlot$abs_seg
  
  rswgsabsolutecn::plotProfile(relcn,ploidy = ploidy,purity = purity,
              clonality = errors["clonality"],rmse = errors["rmse"])
}
dev.off()

