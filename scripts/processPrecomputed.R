args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(QDNAseqmod))
suppressPackageStartupMessages(library(Biobase))

cat("[processPrecomputed] Generating file to skip stage_1\n")
cat("[processPrecomputed] Reading config and samplesheet files...\n")
config <- read_yaml(file="config/config.yaml")

samplesheet <- read.table(config$samplesheet,header=TRUE,sep="\t")
projectBin <- paste0(config$project_name,"_",config$bin,"kb")
outputLoc <- paste0(config$out_dir,"sWGS_fitting/",projectBin,"/")


if(any(is.na(samplesheet$precPloidy)) | any(is.na(samplesheet$precPurity))){
	stop("some samples do not have precomputed ploidy and purity values - check the samplesheet")
}

cat("[processPrecomputed] Loading genome bin data...\n")
bin <- config$bin
bin_size <- bin*1000
bins<-getBinAnnotations(binSize = bin)
# actual number of bins varies in stage_1 due to filtering so total usable bins is adjusted to attempt to fix this
nbins_ref_genome <- round(sum(bins@data$use) * 0.932) #0.932 ratio between actual usable bins and total usable bins in bin annotation data

pre <- "absolute_PRE_down_sampling/"
preFile <- paste0(outputLoc,pre,config$project_name,"_fit_QC_predownsample.tsv")


if(!dir.exists(outputLoc)){
        dir.create(outputLoc,recursive=TRUE)
}

cat("[processPrecomputed] Verifying BAM files...\n")
# check bams
bam <- samplesheet$file
outname <- paste0(outputLoc,"bam.ok")
log_vector <- c()
check_vector <- c()
for(i in bam){
        cmd <- paste0("samtools quickcheck -q ",i)
        if(system(cmd) == 0){
                log_vector <- append(log_vector, paste0("BAM valid - ",i))
                check_vector <- append(check_vector,FALSE)
        } else {
                log_vector <- append(log_vector,paste0("BAM invalid or missing - ",i))
                check_vector <- append(check_vector,TRUE)
        }
}

if(any(check_vector)){
        outname <- gsub(pattern = "ok",replacement = "invalid",outname)
        writeLines(text = as.character(log_vector),con = outname)
} else {
        writeLines(text = as.character(log_vector),con = outname)
}

cat("[processPrecomputed] Creating BAM file symlinks...\n")
# generate symlinks for bams
symoutLoc <- paste0(outputLoc,"bams/")
if(!dir.exists(symoutLoc)){
        dir.create(symoutLoc,recursive=TRUE)
}
for(i in 1:nrow(samplesheet)){
	symcmd <- paste0("ln -s ",samplesheet$file[i]," ",symoutLoc,samplesheet$SAMPLE_ID[i],".bam")
	system(symcmd)
}

cat("[processPrecomputed] Generating spoof Pre-downsampled RDS files...\n")
# spoof RDS files from stage_1
rdsoutLoc <- paste0(outputLoc,pre,"relative_cn_rds/")
if(!dir.exists(rdsoutLoc)){
        dir.create(rdsoutLoc,recursive=TRUE)
}

for(i in 1:nrow(samplesheet)){
	writeLines(text = as.character(samplesheet$SAMPLE_ID[i]),
		paste0(rdsoutLoc,config$project,"_",samplesheet$SAMPLE_ID[i],"_",bin,"kb_relSmoothedCN.rds"))
}
writeLines(text = as.character(samplesheet$SAMPLE_ID),
                paste0(outputLoc,pre,projectBin,"_relSmoothedCN.rds"))

cat("[processPrecomputed] Generating stage_2 QC file input...\n")
# generate qc file
preFileCols <- c("clonality","powered","TP53cn","expected_TP53_AF","TP53freq","rank_clonality","pl_diff","pu_diff","new","new_state")
samplesheet <- samplesheet %>%
		dplyr::select(SAMPLE_ID,PATIENT_ID,precPloidy,precPurity,TP53freq,smooth) %>%
		rename("ploidy"="precPloidy","purity"="precPurity") %>%
		mutate(!!!setNames(rep(NA, length(preFileCols)), preFileCols)) %>%
		mutate(downsample_depth = (((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)) %>%
		mutate(use = TRUE) %>%
		mutate(notes = "using preprocessed ploidy purity values") %>%
		select(SAMPLE_ID,PATIENT_ID,ploidy,purity,clonality,downsample_depth,powered,TP53cn,expected_TP53_AF,TP53freq,smooth,rank_clonality,pl_diff,pu_diff,new,new_state,use,notes)

foutputLoc <- paste0(outputLoc,pre)
if(!dir.exists(foutputLoc)){
	dir.create(foutputLoc,recursive=TRUE)
}

cat("[processPrecomputed] Marking run as precomputed...\n")
write.table(samplesheet,preFile,sep="\t",col.names=T,row.names=F,quote=F)
write.table(samplesheet[,c("SAMPLE_ID","ploidy","purity")],paste0(foutputLoc,"precomp.ok"),
	sep="\t",col.names=T,row.names=F,quote=F)
writeLines(text = c("PRECOMPUTED"),paste0(foutputLoc,"GENERATED_USING_PRECOMPUTED_PL_PU_RDS_FILES_ARE_EMPTY"))

cat("[processPrecomputed] Processing precomputed values complete\n")
cat("[processPrecomputed] stage_2 can now run without completing stage_1\n")
#END
