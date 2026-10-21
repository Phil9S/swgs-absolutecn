# split input
options(scipen = 999)
args = commandArgs(trailingOnly=TRUE)

# Redirect R output, messages, warnings and errors to the Snakemake rule log
if(exists("snakemake")){
  logFile <- file(snakemake@log[[1]], open = "wt")
  sink(logFile, type = "output")
  sink(logFile, type = "message")
}

if(exists("snakemake")){
  tsv <- snakemake@input[["tsv"]]
  tsvOut <- snakemake@output[["tsv"]]
  sampleName <- snakemake@params[["sample"]]
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "splitInput"))
  tsv <- opts$tsv
  tsvOut <- opts$tsvOut
  sampleName <- opts$sampleName
}

QCTable <- read.table(tsv,header = T,sep = "\t",colClasses = c(SAMPLE_ID = "character", PATIENT_ID = "character"))
QCTableSplit <- split(QCTable,f = as.factor(QCTable$SAMPLE_ID))

splitTable <- QCTableSplit[[sampleName]]
write.table(x = splitTable,
	    file = tsvOut,
	    append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,
	    col.names = TRUE)
