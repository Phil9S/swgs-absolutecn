library(yaml)

config <- read_yaml(file="config/config.yaml")

files <- sapply(config$bins,FUN=function(x) paste0(config$out_dir,"sWGS_fitting/",config$project_name,"_",x,"kb/absolute_PRE_down_sampling/",config$project_name,"_fit_QC_predownsample.tsv"))

for(file in files){
  x <- read.table(file = file,header=T,sep="\t")
  x$use <- !duplicated(x$SAMPLE_ID)
  write.table(x,file,row.names=F,col.names=T,sep="\t",quote=F)
}
