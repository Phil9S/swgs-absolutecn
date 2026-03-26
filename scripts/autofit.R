# CN Profile autofit

# grab commandline arguments passed via snakemake object
args = commandArgs(trailingOnly=TRUE)
fitsFile <- snakemake@input[["tsv"]]
fitsTable <- read.table(file = fitsFile,header=T,sep="\t")
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
outputFile <- snakemake@output[["tsv"]]
flagThreshold <- snakemake@params[["flgThreshold"]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]

fitModel <- readRDS("resources/swgs_abs_rfranger_model.rds")

source("scripts/funcs.R")
options(scipen = 999)

fitsTableTriage <- as.data.frame(predictProfile(qctable = fitsTable,
                                                model = fitModel,
                                                method = "randforest",
                                                flagThreshold = 0.84,
                                                multiFitFilter = "clonality"))

write.table(fitsTableTriage,file = outputFile,append = FALSE,
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

