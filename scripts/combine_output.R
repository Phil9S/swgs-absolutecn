## combine_output.R
# Combines sample-level outputs in the pre-downsample phase
# (stage_1) of swgs-abscn to allow for split processing and manual fitting
args = commandArgs(trailingOnly=TRUE)

rds.filename <- snakemake@input[["rds"]]
filelist <- snakemake@input[["tsv"]]

source(file.path(snakemake@scriptdir,"funcs.R"))
options(scipen = 999)

combined.tsv <- do.call(rbind,lapply(filelist,FUN = function(x){
          tab <- read.table(x,sep="\t",header = T)
          return(tab)
}))

rds.list <- lapply(rds.filename,readRDS)
combined.RDS <- collapse_rds(rds.list)

saveRDS(object = combined.RDS,snakemake@output[["rds"]])
write.table(x = combined.tsv,file = snakemake@output[["tsv"]],
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
