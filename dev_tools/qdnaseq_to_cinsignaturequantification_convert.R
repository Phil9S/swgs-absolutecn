## Run CINSignatureQuant on swgs-abscn results
args <- commandArgs(trailingOnly=T)
library(CINSignatureQuantification)
library(tidyverse)
library(yaml)

config <- read_yaml(file="config/config.yaml")

projectBin <- paste0(config$project_name,"_",config$bin,"kb")
outputLoc <- paste0(config$out_dir,"sWGS_fitting/",projectBin,"/")
post <- "absolute_POST_down_sampling/abs_cn_rds/"

cn_file <- paste0(outputLoc,post,projectBin,"_ds_absCopyNumber.rds")
if(file.exists(cn_file)){

  re <- readRDS(cn_file)
  segs <- CINSignatureQuantification:::getSegTable(re)

  segs$start <- as.numeric(segs$start)
  segs$end <- as.numeric(segs$end)
  segs$segVal <- as.numeric(segs$segVal)
  segs$segVal[segs$segVal < 0] <- 0

  q <- CINSignatureQuantification::quantifyCNSignatures(segs,experimentName = "PEO_sc_pseudo_bulk")

  acts <- CINSignatureQuantification::getActivities(q)

  write.table(x = segs,file = paste0(outputLoc,projectBin,"_ds_absCopyNumber_segTable.tsv"),
              sep = "\t",row.names = T,col.names = T,quote = F)
              
  write.table(x = acts,file = paste0(outputLoc,projectBin,"_ds_absCopyNumber_acts.tsv"),
              sep = "\t",row.names = T,col.names = T,quote = F)  
                        
  saveRDS(object = q,file = paste0(outputLoc,projectBin,"_ds_absCopyNumber_cinsigquant_obj.rds"))

} else {
  cat("no downsampled cn profiles found\n")
}

