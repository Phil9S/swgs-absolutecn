# gridsearch_filtering.R 
## Outputs a table of absolute copy number fits across a ploidy and purity gridsearch
## using mean absolute error (aka clonality) error function, sorted by lowest.
## it is designed to be used in conjunction with TP53 allele frequency to determine
## precise purity and ploidy fit
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[["rds"]]
filelist <- snakemake@input[["cl"]]
metafile <- snakemake@params[["meta"]]
metadata <- read.table(file = metafile,header=T,sep="\t")
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]
af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
filter_underpowered <- as.logical(snakemake@params[["filter_underpowered"]]) # Filter for powered fits only
filter_homozygous <- as.logical(snakemake@params[["filter_homozygous"]]) # filter homozygous loss
homozygous_prop <- as.numeric(snakemake@params[["homozygous_prop"]])

# source functions (temp until moved to package) and set scipen
source("scripts/funcs.R")
options(scipen = 999)
outpath <- paste0(out_dir,"sWGS_fitting/",
                  project,"_",bin,"kb/absolute_PRE_down_sampling/")

# read in relative CN data
# collapse rds files function

rds.list <- lapply(rds.filename,FUN=function(x){
  readRDS(x)
})

# Combine and load rds objects - COULD CHANGE TO one-by-one rather than jointly
relative_smoothed <- collapse_rds(rds.list)
saveRDS(relative_smoothed,paste0(outpath,project,"_",bin,"kb_relSmoothedCN.rds"))

filelist <- snakemake@input[["cl"]]
fitTable <- do.call(rbind,
			lapply(filelist,FUN = function(x){
				tab <- read.table(x,sep="\t",skip=1)
				return(tab)
			}))
colnames(fitTable) <- fittingColumnNames

fitTable <- dplyr::left_join(fitTable,metadata,by="SAMPLE_ID") %>%
  dplyr::select(-file) %>%
  dplyr::relocate(PATIENT_ID,.after = SAMPLE_ID) %>%
  dplyr::relocate(TP53freq,smooth,.after = expected_TP53_AF) %>%
  dplyr::relocate(precPloidy,precPurity,.after = purity)

## Apply hard filters
##  filter under powered fits when config variable is TRUE
if(filter_underpowered){
  fitTable <- fitTable %>%
    dplyr::filter(powered == 1) 
}
## filter high prop homozygous loss when config variable is TRUE
if(filter_homozygous){
  fitTable <- fitTable %>%
    dplyr::filter(homozygousLoss <= homozygous_prop)
}

# standard filtering
filtered_results <- fitTable %>%
  #dplyr::select(SAMPLE_ID, PATIENT_ID, everything()) %>% redundant select
  dplyr::group_by(SAMPLE_ID, ploidy) %>%
  dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% #rank clonality within a unique ploidy state 
  dplyr::filter(rank_clonality == 1) %>% #select ploidy with the lowest clonality within a unique ploidy state 
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::top_n(-10, wt = clonality) %>% # select top 10 ploidy states with the lowest clonality values
  dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% # rank by clonality within a sample across ploidy in top 10
  # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
  dplyr::filter(is.na(TP53freq) | dplyr::near(expected_TP53_AF,TP53freq,tol = af_cutoff)) %>%  
  dplyr::arrange(PATIENT_ID, SAMPLE_ID)

# Further limit the results by selecting adjacent ploidy states with the lowest clonality 
# values where multiple similar solutions are present. Plody values within a 
# threshold of 0.3 grouped together and the lowest clonality value is selected
pruned_results <- filtered_results %>%
  dplyr::arrange(SAMPLE_ID, ploidy) %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::mutate(pl_diff = abs(ploidy - dplyr::lag(ploidy))) %>% #, pu_diff = abs(purity - dplyr::lag(purity)) not used
  dplyr::mutate(new_state_n = dplyr::row_number() == 1 | pl_diff > 0.3) %>%
  dplyr::mutate(new_state = cumsum(new_state_n)) %>%
  dplyr::group_by(SAMPLE_ID, new_state) %>%
  dplyr::filter(rank_clonality == min(rank_clonality)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(use = rep(NA,times=nrow(.)),notes = rep(NA,times=nrow(.)))

write.table(filtered_results,paste0(outpath,project,"_filtered_results.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

write.table(pruned_results,paste0(outpath,project,"_fit_QC_predownsample.tsv"),
  sep="\t",col.names=T,row.names=F,quote=F)

## ADDED by PS - adding output folder for results
if(!dir.exists(paste0(outpath,"plots"))){
	dir.create(paste0(outpath,"plots"),recursive = TRUE)
}

for(sample in unique(pruned_results$SAMPLE_ID)){
    dat <- pruned_results %>%
      dplyr::filter(SAMPLE_ID == sample) %>%
      dplyr::arrange(ploidy)
    
    relcn <- relative_smoothed[Biobase::fData(relative_smoothed)$use,sample]
    cn <- Biobase::assayDataElement(relcn,"copynumber")
    seg <- Biobase::assayDataElement(relcn,"segmented")
    
    rel_ploidy <- mean(cn,na.rm=T)
    ll <- nrow(dat)
    
    png(paste0(outpath,"plots/", sample, ".png"),type = "cairo", w= 450*ll, h = 350)
    par(mfrow = c(1,ll)) 
    for(n in 1:nrow(dat)){
      
      ploidy <- dat[n,]$ploidy
      purity <- dat[n,]$purity
      cellploidy <- ploidy * purity + (2*(1-purity))
      seqdepth <- rel_ploidy/cellploidy
      
      abs_cn <- depthtocn(cn,purity,seqdepth)
      abs_seg <- depthtocn(seg,purity,seqdepth)
      
      integer_seg <- round(abs_seg,digits = 0)
      
      errors <- abs_seg - integer_seg
      clonality <- mean(abs(errors)) # clonality is a legacy name for MAE
      rmse <- sqrt(mean(errors^2)) # Root Mean Squared Error
      
      Biobase::assayDataElement(relcn,"copynumber") <- abs_cn
      Biobase::assayDataElement(relcn,"segmented") <- abs_seg
      
      plotProfile(relcn,ploidy = ploidy,purity = purity,
                  clonality = clonality,rmse = rmse)
    }
    dev.off()
}

#END
