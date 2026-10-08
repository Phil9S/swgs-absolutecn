# split input
options(scipen = 999)
args = commandArgs(trailingOnly=TRUE)

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

QCTable <- read.table(tsv,header = T,sep = "\t")
QCTableSplit <- split(QCTable,f = as.factor(QCTable$SAMPLE_ID))

splitTable <- QCTableSplit[[sampleName]]
write.table(x = splitTable,
	    file = tsvOut,
	    append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,
	    col.names = TRUE)
