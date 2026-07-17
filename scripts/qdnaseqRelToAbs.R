## qdnaseq_rel_to_abs.R
# Performs the final conversion of CN profiles from
# relative to absolute fitting space using the selected ploidy/purity combination
# of the CBS-segmented RDS object generated from the downsampled BAM file

args = commandArgs(trailingOnly=TRUE)
options(scipen = 999)
`%>%` <- dplyr::`%>%`

#inputs
rds <- snakemake@input[["rds"]]
meta <- snakemake@input[["meta"]]
metadata <- read.table(meta,header = T,sep = "\t")

# params
bin <- as.numeric(snakemake@params[["bin"]])
project <- snakemake@params[["project"]]
genome <- snakemake@params[["genome"]]
sampleName <- snakemake@params[["sample"]]

#outputs
tsvOut <- snakemake@output[["tsv"]]
rdsOut <- snakemake@output[["rds"]]
tabOut <- snakemake@output[["tab"]]
plotOut <- snakemake@output[["plot"]]

metadata <- metadata[metadata$use == "TRUE",]

# load rds object
rdsRel <- readRDS(rds)
# Add metadata to pheno information
Biobase::pData(rdsRel) <- dplyr::left_join(Biobase::pData(rdsRel),
                                            metadata,by=c("name"="SAMPLE_ID")) %>% as.data.frame()
rownames(Biobase::pData(rdsRel)) <- Biobase::pData(rdsRel)$name

## fix col types
Biobase::pData(rdsRel)$ploidy <- as.numeric(Biobase::pData(rdsRel)$ploidy)
Biobase::pData(rdsRel)$TP53cn <- as.numeric(Biobase::pData(rdsRel)$TP53cn)
Biobase::pData(rdsRel)$purity <- as.numeric(Biobase::pData(rdsRel)$purity)
Biobase::pData(rdsRel)$expected_TP53_AF <- as.numeric(Biobase::pData(rdsRel)$expected_TP53_AF)

## Get target anchor gene segments
gene_bin_seg <- rswgsabsolutecn::getGeneSeg(abs_data = rdsRel,genome = genome)
## Extract CN and Segs 
abs_profiles <- rdsRel[Biobase::fData(rdsRel)$use,]

ploidy <- Biobase::pData(abs_profiles)$ploidy
purity <- Biobase::pData(abs_profiles)$purity
gridVals <- rswgsabsolutecn::gridStats(obj = abs_profiles,
                      ploidy = ploidy,
                      purity = purity)

num_segs <- length(rle(as.numeric(gridVals$seg))$values)
errors <- rswgsabsolutecn::getErrors(data = gridVals)

targetCNVal <- median(gridVals$seg[gene_bin_seg])
TP53cn <- round(rswgsabsolutecn::depthtocn(targetCNVal,purity,gridVals$seqdepth),1) # to 1 decimal place
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
colnames(res) <- get(utils::data("rel2absColumnNames",package = "rswgsabsolutecn",envir = environment()))

dsFittingColumnNames <- get(utils::data("dsFittingColumnNames",package = "rswgsabsolutecn",envir = environment()))
resFormat <- res %>%
  dplyr::select(-c("ploidy","purity")) %>%
  dplyr::left_join(.,Biobase::pData(abs_profiles),
                   by=c("SAMPLE_ID"="name","PATIENT_ID","TP53freq"),
                   suffix = c(".post",".pre")) %>%
  dplyr::relocate(any_of(dsFittingColumnNames)) %>%
  dplyr::select(-c("pl_diff","new_state_n","new_state","rank_clonality")) %>%
  data.frame()

png(plotOut,type="cairo",w = 8,h = 6,unit="in",res = 250)
par(mfrow = c(1,1))
rswgsabsolutecn::plotProfile(abs_profiles,ploidy = ploidy,purity = purity,
            clonality = errors["clonality"],rmse = errors["rmse"])
dev.off()

# Save rds
saveRDS(abs_profiles,file = rdsOut)

# save segTable
segTable <- rswgsabsolutecn::getSegTable(abs_profiles)
write.table(segTable,tabOut,sep = "\t",quote=F,row.names=FALSE)

# write table of fits
write.table(resFormat,tsvOut, sep = "\t",quote=F,row.names=FALSE)

