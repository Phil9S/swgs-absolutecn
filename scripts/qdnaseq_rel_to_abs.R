# QDNAseqmod relative to absolute fitting
args = commandArgs(trailingOnly=TRUE)

rds.filename <- snakemake@input[["rds"]]
metafile <- snakemake@input[["meta"]]
metadata <- read.table(metafile,header = T,sep = "\t")
output_dir <- snakemake@params[["outdir"]]
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
genome <- snakemake@params[["genome"]]
sampleName <- snakemake@params[["sample"]]

source("scripts/funcs.R")

outpath <- paste0(output_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_POST_down_sampling/abs_cn_rds/")
metadata <- metadata[metadata$use == "TRUE",]

## TP53 target bin (genome dependent)
if(genome == "hg19"){
	target <- c("17:7565097-7590863")
} else if(genome == "hg38"){
	target <- c("17:7661779-7687538")
}

if(!dir.exists(paste0(outpath,"plots"))){
	dir.create(paste0(outpath,"plots"))
}

# load rds object
rds.rel <- readRDS(rds.filename)

# Add metadata to pheno information
Biobase::pData(rds.rel) <- dplyr::left_join(Biobase::pData(rds.rel),metadata,by=c("name"="SAMPLE_ID")) %>%
  as.data.frame()
rownames(Biobase::pData(rds.rel)) <- Biobase::pData(rds.rel)$name

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(target = target,abs_data = rds.rel)

# Generate abs plot and table of fits
abs_profiles <- rds.rel[Biobase::fData(rds.rel)$use,]
relcn <- abs_profiles # make mutable copy

# Extract CN and Segs
cn <- Biobase::assayDataElement(relcn,"copynumber")
seg <- Biobase::assayDataElement(relcn,"segmented")
rel_ploidy <- mean(cn,na.rm=T)
ploidy <- Biobase::pData(relcn)$ploidy
purity <- Biobase::pData(relcn)$purity
cellploidy <- ploidy * purity + 2*(1 - purity)
seqdepth <- rel_ploidy / cellploidy
# Convert to abs
abs_cn <- depthtocn(cn,purity,seqdepth)
abs_seg <- depthtocn(seg,purity,seqdepth)

num_segs <- length(rle(abs_seg)$values)
MedianSegVar <- calculateSegmentVar(abs_seg = as.numeric(abs_seg),abs_cn = abs_cn)
integer_seg <- round(abs_seg,digits = 0)
errors <- abs_seg - integer_seg
clonality <- mean(abs(errors)) # clonality is a legacy name for MAE
rmse <- sqrt(mean(errors^2)) # Root Mean Squared Error

Biobase::assayDataElement(relcn,"copynumber") <- abs_cn
Biobase::assayDataElement(relcn,"segmented") <- abs_seg

# Add TP53 info
targetCNVal <- median(seg[gene_bin_seg])
TP53cn <- round(depthtocn(targetCNVal,purity,seqdepth),1) # to 1 decimal place / altered to correct bin value
expected_TP53_AF <- TP53cn * purity / (TP53cn * purity + 2*(1 - purity))
TP53freq <- Biobase::pData(relcn)$TP53freq

# Add patient-level info
pat <- as.character(Biobase::pData(relcn)$PATIENT_ID)
res <- as.data.frame(matrix(c(sampleName,pat,ploidy,purity,num_segs,TP53cn,
                          round(expected_TP53_AF,2),TP53freq,
                          round(clonality,5),
                          round(rmse,5),
                          round(MedianSegVar,5)),nrow = 1,ncol = 11))
# Add to abs RDS
Biobase::assayDataElement(abs_profiles,"copynumber") <- abs_cn
Biobase::assayDataElement(abs_profiles,"segmented") <- abs_seg

png(paste0(outpath,"plots/",sampleName,".png"),type="cairo",w = 8,h = 6,unit="in",res = 250)
par(mfrow = c(1,1))
plotProfile(relcn,ploidy = ploidy,purity = purity,
            clonality = clonality,rmse = rmse)
dev.off()
print("here9")
# Annotated and rename table
colnames(res) <- c("SAMPLE_ID","PATIENT_ID","ploidy","purity","segments","TP53cn",
                   "expected_TP53_AF","TP53freq","clonality","rmse","segvariance")

print("here10")
res
Biobase::pData(abs_profiles)
resFormat <- res %>%
  dplyr::select(-c("ploidy","purity")) %>%
  dplyr::left_join(.,Biobase::pData(abs_profiles),
                   by=c("SAMPLE_ID"="name","PATIENT_ID","TP53freq"),
                   suffix = c(".post",".pre")) %>%
  dplyr::relocate(any_of(dsFittingColumnNames)) %>%
  dplyr::select(-c("pl_diff","new_state_n","new_state"))

resFormat <- data.frame(resFormat,stringsAsFactors = F)
# Save rds
saveRDS(abs_profiles,file=paste0(outpath,project,"_",sampleName,"_",bin,"kb_ds_absCopyNumber.rds"))

# save segTable
segTable <- getSegTable(abs_profiles)

write.table(segTable,paste0(outpath,project,"_",sampleName,"_",bin,"kb_ds_absCopyNumber_segTable.tsv"),
	sep = "\t",quote=F,row.names=FALSE)

#write table of fits
write.table(resFormat,paste0(outpath,project,"_",sampleName,"_",bin,"kb_ds_abs_fits.tsv"),
  sep = "\t",quote=F,row.names=FALSE)
