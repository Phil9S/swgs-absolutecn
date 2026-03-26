# CN Profile autofit

# grab commandline arguments passed via snakemake object
args = commandArgs(trailingOnly=TRUE)
fitsFile <- snakemake@input[["tsv"]]
fitsTable <- read.table(file = fitsFile,header=T,sep="\t")
outputFile <- snakemake@output[["tsv"]]

flagThreshold <- 0.84
fitMethod <- "randforest"
errorMetric <- "clonality"

fitModel <- readRDS("resources/swgs_abs_rfranger_model.rds")

options(scipen = 999)
source("scripts/funcs.R")

fitsTableTriage <- as.data.frame(predictProfile(qctable = fitsTable,
                                                model = fitModel,
                                                method = fitMethod,
                                                flagThreshold = flagThreshold,
                                                errorMetric = errorMetric))

write.table(fitsTableTriage,file = outputFile,append = FALSE,
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

