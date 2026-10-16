## Run CINSignatureQuant on swgs-abscn results
args <- commandArgs(trailingOnly=T)

# Redirect R output, messages, warnings and errors to the Snakemake rule log
if(exists("snakemake")){
  logFile <- file(snakemake@log[[1]], open = "wt")
  sink(logFile, type = "output")
  sink(logFile, type = "message")
}

# Explicit library loading needed for package data as CINSignatureQuantification package
# does not explicitly reference datasets
if (!requireNamespace("CINSignatureQuantification", quietly = TRUE)) {
  stop(
    "\nCINSignatureQuantification is required when computeSigs = TRUE.\n\n",
    "This package is licensed under the GAP Available Source License v1.0 ",
    "(non-commercial academic use only) \nand is NOT bundled with this pipeline.",
    " See https://github.com/markowetzlab/CINSignatureQuantification/blob/main/LICENSE).",
    "\n\n",
    "For virtual environment installations run:\n",
    "Rscript -e 'remotes::install_github(\"markowetzlab/CINSignatureQuantification\")'\n\n",
    "For containerised implementation please see swgs-absolutecn README\n")
}

# Licensed under GAP ASL v1.0 because this script
# depends on and is a derivative interface to CINSignatureQuantification.
library(CINSignatureQuantification)

if(exists("snakemake")){
  tsv <- snakemake@input[["tsv"]]
  tsvOutDrews <- snakemake@output[["tsvDrews"]]
  tsvOutMac <- snakemake@output[["tsvMac"]]
  bin <- as.numeric(snakemake@params[["bin"]])
  project <- snakemake@params[["project"]]
  genome <- as.character(snakemake@params[["genome"]])
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "computeSignatures"))

  tsv <- opts$tsv
  tsvOutDrews <- opts$tsvDrews
  tsvOutMac <- opts$tsvMac
  bin <- as.numeric(opts$bin)
  project <- opts$project
  genome <- opts$genome
}

projectBin <- paste0(project,"_",bin,"kb")

## Load seg table and validate input columns
segs <- read.table(file = tsv,header = TRUE,sep = "\t",
  colClasses = c(sample = "character", chromosome = "character",
  start = "numeric", end = "numeric", segVal = "numeric"))
segs$start <- as.numeric(segs$start)
segs$end <- as.numeric(segs$end)
segs$segVal <- as.numeric(segs$segVal)
# not required due to pipeline validation, but just in case it is 
segs$segVal[segs$segVal < 0] <- 0
# Generate signature activities using Drews and Mac methods
cnsigs_drews <- CINSignatureQuantification::quantifyCNSignatures(segs,
  experimentName = projectBin,
  method = "drews",build = genome)
cnsigs_drews_acts <- CINSignatureQuantification::getActivities(cnsigs_drews)

cnsigs_mac <- CINSignatureQuantification::quantifyCNSignatures(segs,
  experimentName = projectBin,
  method = "mac",build = genome)
cnsigs_mac_acts <- CINSignatureQuantification::getActivities(cnsigs_mac)

write.table(x = cnsigs_drews_acts,file = tsvOutDrews,
  sep = "\t",row.names = T,col.names = T,quote = F)

saveRDS(object = cnsigs_drews,file =paste0(dirname(tsvOutDrews),"/",
  gsub(x = basename(tsvOutDrews),pattern = "_acts.tsv",replacement = "_obj.rds")))

write.table(x = cnsigs_mac_acts,file = tsvOutMac,
  sep = "\t",row.names = T,col.names = T,quote = F)

saveRDS(object = cnsigs_mac,file =paste0(dirname(tsvOutMac),"/",
  gsub(x = basename(tsvOutMac),pattern = "_acts.tsv",replacement = "_obj.rds")))
