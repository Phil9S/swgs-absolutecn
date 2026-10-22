## combine_output.R
# Combines sample-level outputs (stage_1 and stage_2)
# of swgs-abscn to match previous (<2.0.0) version using staged implementation
args = commandArgs(trailingOnly=TRUE)

# Redirect R output, messages, warnings and errors to the Snakemake rule log
if(exists("snakemake")){
  logFile <- file(snakemake@log[[1]], open = "wt")
  sink(logFile, type = "output")
  sink(logFile, type = "message")
}

if(exists("snakemake")){
  tsv <- snakemake@input[["tsv"]]
  rds <- snakemake@input[["rds"]]
  tab <- snakemake@input[["tab"]]
  tsvOut <- snakemake@output[["tsv"]]
  rdsOut <- snakemake@output[["rds"]]
  tabOut <- snakemake@output[["tab"]]
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "combineResults"))
  
  # comma-separated lists of per-sample files
  tsv <- strsplit(opts$tsv, ",")[[1]]
  rds <- strsplit(opts$rds, ",")[[1]]
  tab <- if(is.null(opts$tab)) NULL else strsplit(opts$tab, ",")[[1]]
  tsvOut <- opts$tsvOut
  rdsOut <- opts$rdsOut
  tabOut <- opts$tabOut
}

options(scipen = 999)

## combining stage_2 outputs
if(!is.null(tab)){
  combinedTab <- do.call(rbind,lapply(tab,FUN = function(x){
  y <- read.table(x,sep="\t",header = T,
  colClasses = c(sample = "character"))
  return(y)
}))

write.table(x = combinedTab,file = tabOut,
  append = F,quote = F,
  sep = "\t",row.names = F,col.names = T)
}


combinedTsv <- do.call(rbind,lapply(tsv,FUN = function(x){
  y <- read.table(x,sep="\t",header = T,
  colClasses = c(SAMPLE_ID = "character", PATIENT_ID = "character"))
  return(y)
}))

rdsList <- lapply(rds,readRDS)
combinedRDS <- rswgsabsolutecn::collapseRDS(rdsList)

saveRDS(object = combinedRDS,rdsOut)

write.table(x = combinedTsv,file = tsvOut,
  append = F,quote = F,sep = "\t",
  row.names = F,col.names = T)
