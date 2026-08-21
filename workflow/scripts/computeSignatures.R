## Run CINSignatureQuant on swgs-abscn results
args <- commandArgs(trailingOnly=T)

# Explicit library loading needed for package data as CINSignatureQuantification package
# does not explicitly reference datasets
library(CINSignatureQuantification)

if(exists("snakemake")){
  tsv <- snakemake@input[["tsv"]]
  tsvOutDrews <- snakemake@output[["tsvDrews"]]
  tsvOutMac <- snakemake@output[["tsvMac"]]
  bin <- as.numeric(snakemake@params[["bin"]])
  project <- snakemake@params[["project"]]
  genome <- as.character(snakemake@params[["genome"]])
  sampleName <- as.character(snakemake@params[["sample"]])
} else {
  opts <- optparse::parse_args(rswgsabsolutecn::setArgs(script = "computeSignatures"))

  tsv <- opts$tsv
  tsvOutDrews <- opts$tsvDrews
  tsvOutMac <- opts$tsvMac
  bin <- as.numeric(opts$bin)
  project <- opts$project
  genome <- opts$genome
  sampleName <- as.character(opts$sampleName)
}

projectBin <- paste0(project,"_",bin,"kb")

## Load seg table and validate input columns
segs <- read.table(file = tsv,header = TRUE,sep = "\t",as.is = TRUE)
segs$start <- as.numeric(segs$start)
segs$end <- as.numeric(segs$end)
segs$segVal <- as.numeric(segs$segVal)
segs$segVal[segs$segVal < 0] <- 0
## Generate
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
                                         gsub(x = basename(tsvOutDrews),
                                              pattern = "_acts.tsv",replacement = "_obj.rds")))

write.table(x = cnsigs_mac_acts,file = tsvOutMac,
            sep = "\t",row.names = T,col.names = T,quote = F)

saveRDS(object = cnsigs_mac,file =paste0(dirname(tsvOutMac),"/",
                                        gsub(x = basename(tsvOutMac),
                                             pattern = "_acts.tsv",replacement = "_obj.rds")))
