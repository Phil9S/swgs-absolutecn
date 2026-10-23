# processPrecomputed.R
args <- commandArgs(trailingOnly=TRUE)

`%>%` <- dplyr::`%>%`

cat("[processPrecomputed] Generating file to skip stage_1\n")
cat("[processPrecomputed] Reading config and samplesheet files...\n")
config <- yaml::read_yaml(file="config/config.yaml")
samplesheet <- read.table(config$samplesheet,header=TRUE,sep="\t",quote = "",comment.char = "",
  colClasses = c(SAMPLE_ID = "character", PATIENT_ID = "character"))
projectBin <- paste0(config$project_name,"_",config$bins,"kb")
# Ensure single trailing slash on out_dir
outDir <- sub("/*$","/",config$out_dir)
outputLoc <- paste0(outDir,"sWGS_fitting/",projectBin,"/")
genome <- as.character(config$genome)

if(any(is.na(samplesheet$precPloidy)) | any(is.na(samplesheet$precPurity))){
  stop(paste0("Missing precomputed ploidy and purity values - check ",config$samplesheet))
}

cat("[processPrecomputed] Loading genome bin data...\n")
bin <- config$bins
bin_size <- bin * 1000
bins <- QDNAseqmod::getBinAnnotations(binSize = bin, genome = genome)

# Actual number of bins varies in stage_1 due to filtering so total usable bins
# is adjusted to attempt to fix this. 0.932 ratio between actual usable bins and
# total usable bins in bin annotation data
nbins_ref_genome <- round(sum(bins@data$use) * 0.932)
project <- config$project_name
pre <- "absolute_PRE_down_sampling/"
preFile <- paste0(outputLoc,pre,project,"_",bin,"kb_combined_QC_predownsample.tsv")

if(!dir.exists(outputLoc)){
  dir.create(outputLoc,recursive=TRUE)
}

symoutLoc <- paste0(outputLoc,"bams/")
if(!dir.exists(symoutLoc)){
  dir.create(symoutLoc,recursive=TRUE)
}

cat("[processPrecomputed] Verifying BAM files...\n")
missingBams <- samplesheet$SAMPLE_ID[!file.exists(samplesheet$file)]
if(length(missingBams) > 0){
  stop(paste0("BAM file(s) not found for sample(s): ",paste0(missingBams,collapse = ", "),
              " - check the 'file' column in ",config$samplesheet))
}

okFiles <- paste0(symoutLoc,samplesheet$SAMPLE_ID,"_bam.ok")
for(i in 1:nrow(samplesheet)){
  rswgsabsolutecn::bamCheck(x = samplesheet$file[i],filetype = config$filetype,
                            outname = okFiles[i])
}
failedBams <- samplesheet$SAMPLE_ID[!file.exists(okFiles)]
if(length(failedBams) > 0){
  stop(paste0("BAM check failed for sample(s): ",paste0(failedBams,collapse = ", "),
              " - see *_bam.invalid files in ",symoutLoc))
}

cat("[processPrecomputed] Creating BAM file symlinks...\n")
# generate symlinks for BAMs
for(i in 1:nrow(samplesheet)){
  bamAbs <- normalizePath(samplesheet$file[i], mustWork = TRUE)
  linkOut <- paste0(symoutLoc,samplesheet$SAMPLE_ID[i],".bam")
  if(file.exists(linkOut) || nzchar(Sys.readlink(linkOut))){
    unlink(linkOut)
  }
  if(!file.symlink(bamAbs, linkOut)){
    stop(paste0("Failed to create symlink: ",linkOut," -> ",bamAbs))
  }
}


cat("[processPrecomputed] Generating spoof Pre-downsampled RDS files...\n")
# spoof RDS files from stage_1
rdsoutLoc <- paste0(outputLoc,pre,"relative_cn_rds/")
if(!dir.exists(rdsoutLoc)){
  dir.create(rdsoutLoc,recursive=TRUE)
}

for(i in 1:nrow(samplesheet)){
  writeLines(text = as.character(samplesheet$SAMPLE_ID[i]),
    paste0(rdsoutLoc,project,"_",samplesheet$SAMPLE_ID[i],"_",bin,"kb_relSmoothedCN.rds"))
}

writeLines(text = as.character(samplesheet$SAMPLE_ID),
  paste0(outputLoc,pre,"relative_cn_rds/",projectBin,"_relSmoothedCN.rds"))

# generate QC file
cat("[processPrecomputed] Generating stage_2 QC file input...\n")
fittingColumnNames <- get(utils::data("fittingColumnNames",package = "rswgsabsolutecn",envir = environment()))
preFileCols <- c(fittingColumnNames,"PATIENT_ID","precPloidy","precPurity","TP53freq","smooth",
"rank_clonality","pl_diff","new_state_n","new_state")

exclude <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","smooth","precPloidy","precPurity","TP53freq")
samplesheetNew <- samplesheet %>%
  dplyr::mutate(ploidy = precPloidy,purity = precPurity) %>%
  dplyr::relocate(ploidy,purity,.before="precPloidy") %>%
  dplyr::mutate(!!!setNames(rep(NA, length(preFileCols)), preFileCols[which(!preFileCols %in% exclude)])) %>%
  dplyr::mutate(downsample_depth = rswgsabsolutecn::getDownsampleDepth(ploidy=ploidy,purity=purity,nbins_ref_genome=nbins_ref_genome)) %>%
  dplyr::mutate(use = TRUE) %>%
  dplyr::mutate(notes = "using preprocessed ploidy purity values") %>%
  dplyr::select(all_of(preFileCols),use,notes)

foutputLoc <- paste0(outputLoc,pre)
if(!dir.exists(foutputLoc)){
  dir.create(foutputLoc,recursive=TRUE)
}

cat("[processPrecomputed] Marking run as precomputed...\n")
write.table(samplesheetNew,preFile,sep="\t",col.names=T,row.names=F,quote=F)
write.table(samplesheetNew[,c("SAMPLE_ID","ploidy","purity")],paste0(foutputLoc,"precomp.ok"),
  sep="\t",col.names=T,row.names=F,quote=F)
writeLines(text = c("PRECOMPUTED"),paste0(foutputLoc,"GENERATED_USING_PRECOMPUTED_PL_PU_RDS_FILES_ARE_EMPTY"))

cat("[processPrecomputed] Processing precomputed values complete\n")
cat("[processPrecomputed] stage_2 can now run without completing stage_1\n")
#END
