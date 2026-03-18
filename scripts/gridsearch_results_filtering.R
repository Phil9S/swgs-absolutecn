# gridsearch_filtering.R 
## Outputs a table of absolute copy number fits across a ploidy and purity gridsearch
## using mean absolute error (aka clonality) error function, sorted by lowest.
## it is designed to be used in conjunction with TP53 allele frequency to determine
## precise purity and ploidy fit
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[["rds"]]
filelist <- snakemake@input[["cl"]]
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
sampleName <- as.character(snakemake@params[["sample"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
filter_underpowered <- as.logical(snakemake@params[["filter_underpowered"]]) # Filter for powered fits only
filter_homozygous <- as.logical(snakemake@params[["filter_homozygous"]]) # filter homozygous loss
homozygous_prop <- as.numeric(snakemake@params[["homozygous_prop"]])

# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)
outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_PRE_down_sampling/")

# read in relative CN data
relative_smoothed <- readRDS(rds.filename)

fitTable <- read.table(filelist,sep="\t",header = TRUE)
colnames(fitTable) <- fittingColumnNames

filteredTables <- filterFitTable(table = fitTable,metadata = metadata,
                                 filter_underpowered = filter_underpowered,
                                 filter_homozygous = filter_homozygous,
                                 af_cutoff = af_cutoff)

write.table(filteredTables$filtered,paste0(outpath,project,"_",sampleName,"_filtered_results.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

write.table(filteredTables$pruned,paste0(outpath,project,"_",sampleName,"_fit_QC_predownsample.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

## ADDED by PS - adding output folder for results
if(!dir.exists(paste0(outpath,"plots"))){
	dir.create(paste0(outpath,"plots"),recursive = TRUE)
}

#for(sample in unique(pruned_results$SAMPLE_ID)){
#dat <- pruned_results %>%
  #dplyr::filter(SAMPLE_ID == unique(pruned_results$SAMPLE_ID)) %>%
  #dplyr::arrange(ploidy)

to_use <- Biobase::fData(relative_smoothed)$use
relcn <- relative_smoothed[to_use,]
cn <- Biobase::assayDataElement(relcn,"copynumber")
seg <- Biobase::assayDataElement(relcn,"segmented")

rel_ploidy <- mean(cn,na.rm=T)
ll <- nrow(pruned_results)

png(paste0(outpath,"plots/",sampleName,".png"),type = "cairo", w= 450*ll, h = 350)
par(mfrow = c(1,ll)) 
for(n in 1:nrow(pruned_results)){
  
  ploidy <- pruned_results[n,]$ploidy
  purity <- pruned_results[n,]$purity
  cellploidy <- ploidy * purity + (2*(1-purity))
  seqdepth <- rel_ploidy/cellploidy
  
  abs_cn <- depthtocn(cn,purity,seqdepth)
  abs_seg <- depthtocn(seg,purity,seqdepth)
  
  integer_seg <- round(abs_seg,digits = 0)
  
  errors <- abs_seg - integer_seg
  clonality <- mean(abs(errors)) # clonality is a legacy name for MAE
  rmse <- sqrt(mean(errors^2)) # Root Mean Squared Error
  
  Biobase::assayDataElement(relcn,"copynumber") <- abs_cn
  Biobase::assayDataElement(relcn,"segmented") <- abs_seg
  
  plotProfile(relcn,ploidy = ploidy,purity = purity,
              clonality = clonality,rmse = rmse)
}
dev.off()

#END
