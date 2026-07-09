## Dev script to summarise segment counts in both pre and post downsampling

args <- commandArgs(trailingOnly=T)

countSegs <- function(x){
  suppressMessages(library(QDNAseqmod))
  suppressMessages(library(Biobase))

  filer <- readRDS(x)
  pS <- filer[Biobase::featureData(filer)$use]
  pSegs <- apply(Biobase::assayDataElement(pS,"segmented"),
                 MARGIN=2,function(x) length(rle(x)$lengths))
  return(pSegs)
}

cat("report segments - use 'report_seg_counts.R all' for individual seg counts\n")

config <- yaml::read_yaml(file="config/config.yaml")

projectBin <- paste0(config$project_name,"_",config$bin,"kb")
outputLoc <- paste0(config$out_dir,"sWGS_fitting/",projectBin,"/")

pre <- "absolute_PRE_down_sampling"
post <- "absolute_POST_down_sampling"

preFiles <- list.files(path = paste0(outputLoc,pre),
                       recursive = T,full.names = T,
                       pattern = "_relSmoothedCN.rds")

postFiles <- list.files(path = paste0(outputLoc,post),
                       recursive = T,full.names = T ,
                       pattern = "_ds_absCopyNumber.rds")

verbose <- FALSE
if(length(args) > 0){
  if(args[1] == "all"){
    verbose <- TRUE
  }
}

if(any(file.exists(postFiles))){
  cat("\ncounting pre-downsampled segments")
  preSegs <- unlist(lapply(preFiles,FUN = countSegs))
  if(verbose){
    print(preSegs)
  } else {
    cat("\npre-downsampled segments")
    print(summary(preSegs))
  }
} else {
  cat("no pre or post downsampled files found\n")
}

if(any(file.exists(postFiles))){
  cat("\ncounting post-downsampled segments")
  postSegs <- unlist(lapply(postFiles[postFiles],countSegs))
  if(verbose){

    print(postSegs)
  } else {
    cat("\npost-downsampled segments")
    print(summary(postSegs))
  }
} else {
  cat("\nNo post downsampled files found\n")
}
