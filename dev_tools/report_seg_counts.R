args <- commandArgs(trailingOnly=T)
library(yaml)

cat("report segments - use 'report_seg_counts.R all' for individual seg counts\n")

config <- read_yaml(file="config/config.yaml")

projectBin <- paste0(config$project_name,"_",config$bin,"kb")
outputLoc <- paste0(config$out_dir,"sWGS_fitting/",projectBin,"/")

pre <- "absolute_PRE_down_sampling/"
post <- "absolute_POST_down_sampling/abs_cn_rds/"

preFile <- paste0(outputLoc,pre,projectBin,"_relSmoothedCN.rds")
postFile <- paste0(outputLoc,post,projectBin,"_ds_absCopyNumber.rds")

verbose <- FALSE
if(length(args) > 0){
  if(args[1] == "all"){
    verbose <- TRUE
  }
}

if(file.exists(preFile)){
  suppressMessages(library(QDNAseqmod))
  suppressMessages(library(Biobase))
  preS <- readRDS(preFile)
  preS <- preS[featureData(preS)$use]
  preSegs <- apply(assayDataElement(preS,"segmented"),MARGIN=2,function(x) length(rle(x)$lengths))
  cat("\nPre-downsampled segments\n")
  if(verbose){
    print(preSegs)
  } else {
    print(summary(preSegs))
  }
  if(file.exists(postFile)){
    postS <- readRDS(postFile)
    postS <- postS[featureData(postS)$use]
    postSegs <- apply(assayDataElement(postS,"segmented"),MARGIN=2,function(x) length(rle(x)$lengths))
    cat("\nPost-downsampled segments\n")
    if(verbose){
      print(postSegs)
    } else {
      print(summary(postSegs))
    }
  }
} else {
  cat("no pre or post downsampled files found\n")
}
