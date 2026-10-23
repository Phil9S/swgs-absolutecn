library(yaml)

config <- read_yaml(file="config/config.yaml")

project <- config$project_name
# Ensure single trailing slash on out_dir
outDir <- sub("/*$","/",config$out_dir)

files <- sapply(config$bins,FUN=function(x) paste0(outDir,"sWGS_fitting/",project,"_",x,"kb/absolute_PRE_down_sampling/",project,"_",x,"kb_combined_QC_predownsample.tsv"))

for(file in files){
  x <- read.table(file = file,header=T,sep="\t",quote = "",comment.char = "",
    colClasses = c(SAMPLE_ID = "character", PATIENT_ID = "character"))
  x$use <- !duplicated(x$SAMPLE_ID)
  write.table(x,file,row.names=F,col.names=T,sep="\t",quote=F)
}
