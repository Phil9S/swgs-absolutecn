## autofit.R
# automates CN profile fitting using either randomforest or
# min error selection of best absolute copy number profile fit. Additionally
# flags fits which are low confidence or potentially poor
args = commandArgs(trailingOnly=TRUE)

if(exists("snakemake")){
  tsv <- snakemake@input[["tsv"]]
  tsvOut <- snakemake@output[["tsv"]]
  flagThreshold <- snakemake@params[["flagThreshold"]]
  fitMethod <- snakemake@params[["fitMethod"]]
  errorMetric <- snakemake@params[["errorMetric"]]
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "autofit"))
  
  tsv <- opts$tsv
  tsvOut <- opts$tsvOut
  flagThreshold <- opts$flagThreshold
  fitMethod <- opts$fitMethod
  errorMetric <- opts$errorMetric
}

fitsTable <- read.table(file = tsv,header=T,sep="\t")

if(fitMethod == "randforest"){
  fitModel <- get(utils::data("swgs_abs_rfranger_model",
    package = "rswgsabsolutecn",envir = environment()))
} else {
  fitModel <- NULL  
}

options(scipen = 999)

fitsTableTriage <- as.data.frame(rswgsabsolutecn::predictProfile(qctable = fitsTable,
                                                model = fitModel,
                                                method = fitMethod,
                                                flagThreshold = flagThreshold,
                                                errorMetric = errorMetric))

write.table(fitsTableTriage,file = tsvOut,append = FALSE,
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
