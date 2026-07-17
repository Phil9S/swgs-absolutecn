# split input
options(scipen = 999)
args = commandArgs(trailingOnly=TRUE)

tsv <- snakemake@input[["tsv"]]
tsvOut <- snakemake@output[["tsv"]]
sampleName <- snakemake@params[["sample"]]

QCTable <- read.table(tsv,header = T,sep = "\t")
QCTableSplit <- split(QCTable,f = as.factor(QCTable$SAMPLE_ID))

splitTable <- QCTableSplit[[sampleName]]
write.table(x = splitTable,
	    file = tsvOut,
	    append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,
	    col.names = TRUE)
