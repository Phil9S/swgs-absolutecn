# split input
args = commandArgs(trailingOnly=TRUE)

file <- snakemake@input[["tsv"]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
sampleName <- snakemake@params[["sample"]]
outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_POST_down_sampling/")

source("scripts/funcs.R")
options(scipen = 999)

QCTable <- read.table(file,header = T,sep = "\t")
QCTableSplit <- split(QCTable,f = as.factor(QCTable$SAMPLE_ID))

splitTable <- QCTableSplit[[sampleName]]
write.table(x = splitTable,
	    file = paste0(outpath,project,"_",sampleName,"_fitted_QC.tsv"),
	    append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,
	    col.names = TRUE)

#END
