# combine single outputs to one file
args = commandArgs(trailingOnly=TRUE)

rds.filename <- snakemake@input[["rds"]]
filelist <- snakemake@input[["tsv"]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]

outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_PRE_down_sampling/")

source("scripts/funcs.R")
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