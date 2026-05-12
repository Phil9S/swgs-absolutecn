## autofit.R
# automates CN profile fitting using either randomforest or
# min error selection of best absolute copy number profile fit. Additionally
# flags fits which are low confidence or potentially poor
args = commandArgs(trailingOnly=TRUE)
fitsFile <- snakemake@input[["tsv"]]
fitsTable <- read.table(file = fitsFile,header=T,sep="\t")
outputFile <- snakemake@output[["tsv"]]

flagThreshold <- snakemake@params[["flagThreshold"]]
fitMethod <- snakemake@params[["fitMethod"]]
errorMetric <- snakemake@params[["errorMetric"]]

if(fitMethod == "randforest"){
  fitModel <- readRDS("resources/swgs_abs_rfranger_model.rds")
  
} else {
  fitModel <- NULL  
}

options(scipen = 999)
source(file.path(snakemake@scriptdir,"funcs.R"))

fitsTableTriage <- as.data.frame(predictProfile(qctable = fitsTable,
                                                model = fitModel,
                                                method = fitMethod,
                                                flagThreshold = flagThreshold,
                                                errorMetric = errorMetric))

write.table(fitsTableTriage,file = outputFile,append = FALSE,
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
