# split input
args = commandArgs(trailingOnly=TRUE)

file <- snakemake@input[["tsv"]]
outfile <- snakemake@output[["tsv"]]
sampleName <- snakemake@params[["sample"]]

source(file.path(snakemake@scriptdir,"funcs.R"))
options(scipen = 999)

QCTable <- read.table(file,header = T,sep = "\t")
QCTableSplit <- split(QCTable,f = as.factor(QCTable$SAMPLE_ID))

splitTable <- QCTableSplit[[sampleName]]
write.table(x = splitTable,
	    file = outfile,
	    append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,
	    col.names = TRUE)
