## qdnaseq_rel_to_abs.R
# Performs the final conversion of CN profiles from
# relative to absolute fitting space using the selected ploidy/purity combination
# of the CBS-segmented RDS object generated from the downsampled BAM file
args = commandArgs(trailingOnly=TRUE)

#inputs
rds.filename <- snakemake@input[["rds"]]
metafile <- snakemake@input[["meta"]]
metadata <- read.table(metafile,header = T,sep = "\t")

# params
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
genome <- snakemake@params[["genome"]]
sampleName <- snakemake@params[["sample"]]

#outputs
tsv_output <- snakemake@output[["tsv"]]
rds_output <- snakemake@output[["rds"]]
seg_output <- snakemake@output[["seg"]]
plot_output <- snakemake@output[["plot"]]

source(file.path(snakemake@scriptdir,"funcs.R"))
options(scipen = 999)

metadata <- metadata[metadata$use == "TRUE",]

# load rds object
rds.rel <- readRDS(rds.filename)
# Add metadata to pheno information
Biobase::pData(rds.rel) <- dplyr::left_join(Biobase::pData(rds.rel),
                                            metadata,by=c("name"="SAMPLE_ID")) %>%
  as.data.frame()
rownames(Biobase::pData(rds.rel)) <- Biobase::pData(rds.rel)$name

## fix col types
Biobase::pData(rds.rel)$ploidy <- as.numeric(Biobase::pData(rds.rel)$ploidy)
Biobase::pData(rds.rel)$TP53cn <- as.numeric(Biobase::pData(rds.rel)$TP53cn)
Biobase::pData(rds.rel)$purity <- as.numeric(Biobase::pData(rds.rel)$purity)
Biobase::pData(rds.rel)$expected_TP53_AF <- as.numeric(Biobase::pData(rds.rel)$expected_TP53_AF)

#Get target anchor gene segments
gene_bin_seg <- get_gene_seg(abs_data = rds.rel,
                             genome = genome)
# # Extract CN and Segs 
abs_profiles <- rds.rel[Biobase::fData(rds.rel)$use,]

ploidy <- Biobase::pData(abs_profiles)$ploidy
purity <- Biobase::pData(abs_profiles)$purity
gridVals <- gridStats(obj = abs_profiles,
                      ploidy = ploidy,
                      purity = purity)

num_segs <- length(rle(as.numeric(gridVals$seg))$values)
errors <- getErrors(data = gridVals)

targetCNVal <- median(gridVals$seg[gene_bin_seg])
TP53cn <- round(depthtocn(targetCNVal,purity,gridVals$seqdepth),1) # to 1 decimal place
expected_TP53_AF <- TP53cn * purity / (TP53cn * purity + 2 * (1-purity))
TP53freq <- Biobase::pData(abs_profiles)$TP53freq

Biobase::assayDataElement(abs_profiles,"copynumber") <- gridVals$abs_cn
Biobase::assayDataElement(abs_profiles,"segmented") <- gridVals$abs_seg

# Add patient-level info
pat <- as.character(Biobase::pData(abs_profiles)$PATIENT_ID)
res <- as.data.frame(matrix(c(sampleName,pat,ploidy,purity,num_segs,TP53cn,
                          round(expected_TP53_AF,2),TP53freq,
                          round(errors,5)),nrow = 1,ncol = 11))
# Annotated and rename table
colnames(res) <- rel2absColumnNames

resFormat <- res %>%
  dplyr::select(-c("ploidy","purity")) %>%
  dplyr::left_join(.,Biobase::pData(abs_profiles),
                   by=c("SAMPLE_ID"="name","PATIENT_ID","TP53freq"),
                   suffix = c(".post",".pre")) %>%
  dplyr::relocate(any_of(dsFittingColumnNames)) %>%
  dplyr::select(-c("pl_diff","new_state_n","new_state","rank_clonality")) %>%
  data.frame()

png(plot_output,type="cairo",w = 8,h = 6,unit="in",res = 250)
par(mfrow = c(1,1))
plotProfile(abs_profiles,ploidy = ploidy,purity = purity,
            clonality = errors["clonality"],rmse = errors["rmse"])
dev.off()

# Save rds
saveRDS(abs_profiles,file = rds_output)

# save segTable
segTable <- getSegTable(abs_profiles)
write.table(segTable,seg_output,sep = "\t",quote=F,row.names=FALSE)

# write table of fits
write.table(resFormat,tsv_output, sep = "\t",quote=F,row.names=FALSE)
