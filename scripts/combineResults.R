## combine_output.R
# Combines sample-level outputs (stage_1 and stage_2)
# of swgs-abscn to match previous (<2.0.0) version using staged implementation
args = commandArgs(trailingOnly=TRUE)

options(scipen = 999)

rds <- snakemake@input[["rds"]]
tsv <- snakemake@input[["tsv"]]
tab <- snakemake@input[["tab"]]

rdsOut <- snakemake@output[["rds"]]
tsvOut <- snakemake@output[["tsv"]]
tabOut <- snakemake@output[["tab"]]

## combining stage_2 outputs
if(!is.null(tab)){
    combinedTab <- do.call(rbind,lapply(tab,FUN = function(x){
	y <- read.table(x,sep="\t",header = T)
	return(y)
    }))

    write.table(x = combinedTab,file = tabOut,
		            append = F,quote = F,
			    sep = "\t",row.names = F,col.names = T)

}


combinedTsv <- do.call(rbind,lapply(tsv,FUN = function(x){
	y <- read.table(x,sep="\t",header = T)
	return(y)
}))

rdsList <- lapply(rds,readRDS)
combinedRDS <- rswgsabsolutecn::collapseRDS(rdsList)

saveRDS(object = combinedRDS,rdsOut)

write.table(x = combinedTsv,file = tsvOut,
            append = F,quote = F,sep = "\t",
	    row.names = F,col.names = T)
