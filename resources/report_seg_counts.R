## Dev script to summarise segment counts in both pre and post downsampling

args <- commandArgs(trailingOnly=T)

countSegsTab <- function(x,stg="pre"){
  filer <- read.table(x,nrows = 1,header = T,sep = "\t",quote = "",comment.char = "",
    colClasses = c(SAMPLE_ID = "character"))
  pSample <- filer$SAMPLE_ID
  pSegs <- ifelse(stg=="pre",filer$segments,filer$segments.post)
  names(pSegs) <- pSample
  return(pSegs)
}

cat("report_seg_counts.R - use 'report_seg_counts.R all' for individual seg counts\n")

config <- yaml::read_yaml(file="config/config.yaml")

projectBin <- paste0(config$project_name,"_",config$bins,"kb")
# Ensure single trailing slash on out_dir
outDir <- sub("/*$","/",config$out_dir)
outputLoc <- paste0(outDir,"sWGS_fitting/",projectBin,"/")

pre <- "absolute_PRE_down_sampling"
post <- "absolute_POST_down_sampling"

preFilesTab <- list.files(path = paste0(outputLoc,pre),
                       recursive = T,full.names = T,
                       pattern = "_fit_QC_predownsample\\.tsv$")

postFilesTab <- list.files(path = paste0(outputLoc,post),
                        recursive = T,full.names = T ,
                        pattern = "_ds_abs_fits\\.tsv$")
# per-sample tables only - drop the stage_2 combined table
postFilesTab <- postFilesTab[basename(postFilesTab) != paste0(projectBin,"_ds_abs_fits.tsv")]

verbose <- FALSE

if(length(args) > 0){
  if(args[1] == "all"){
    verbose <- TRUE
  }
}

if(any(file.exists(preFilesTab))){
  preSegsTab <- unlist(lapply(preFilesTab,FUN = countSegsTab,stg="pre"))
  if(verbose){
    cat("\npre-downsampled segments\n")
    print(preSegsTab)
  } else {
    cat("\npre-downsampled segments\n")
    print(summary(preSegsTab))
  }
} else {
  cat("no pre downsampled files found\n")
}

if(any(file.exists(postFilesTab))){
  postSegs <- unlist(lapply(postFilesTab,countSegsTab,stg="post"))
  if(verbose){
    cat("\npost-downsampled segments\n")
    print(postSegs)
  } else {
    cat("\npost-downsampled segments\n")
    print(summary(postSegs))
  }
} else {
  cat("\nNo post downsampled files found\n")
}
