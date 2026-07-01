## combine_output.R
# Combines sample-level outputs (stage_1 and stage_2)
# of swgs-abscn to match previous (<2.0.0) version using staged implementation
args = commandArgs(trailingOnly=TRUE)

source(file.path(snakemake@scriptdir,"funcs.R"))
options(scipen = 999)

rds.filename <- snakemake@input[["rds"]]
filelist <- snakemake@input[["tsv"]]

## Only run this if combining stage_2 outputs
if(!is.null(snakemake@input[["tab"]])){
    tablist <- snakemake@input[["tab"]]
    combined.tab <- do.call(rbind,lapply(tablist,FUN = function(x){
	tab <- read.table(x,sep="\t",header = T)
	return(tab)
    }))

    write.table(x = combined.tab,file = snakemake@output[["tab"]],
		            append = F,quote = F,sep = "\t",row.names = F,col.names = T)

}


combined.tsv <- do.call(rbind,lapply(filelist,FUN = function(x){
	tab <- read.table(x,sep="\t",header = T)
	return(tab)
}))

rds.list <- lapply(rds.filename,readRDS)
combined.RDS <- collapse_rds(rds.list)

saveRDS(object = combined.RDS,snakemake@output[["rds"]])
write.table(x = combined.tsv,file = snakemake@output[["tsv"]],
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
