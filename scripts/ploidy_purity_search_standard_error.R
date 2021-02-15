#this script outputs a ranked list of absolute copy number fits using the 
#traditional clonality error function. it is designed to be used in conjuction 
#with TP53 allele frequency to determine precise purity and ploidy fit

#grab commandline arguments
args = commandArgs(trailingOnly=TRUE)
rds.filename <- snakemake@input[[1]]
bin <- as.numeric(snakemake@params[["bin"]])
out_dir <- snakemake@params[["outdir"]]
project <- snakemake@params[["project"]]

#load libraries
suppressPackageStartupMessages(library(QDNAseqmod))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(matrixStats))
suppressWarnings(library(doMC))
suppressWarnings(library(foreach))
registerDoMC(cores)
options(scipen = 999)

#read in relative copy number
rds.obj <- readRDS(rds.filename)
sampleNames(rds.obj)<-str_replace(sampleNames(rds.obj),"[(]","_")
sampleNames(rds.obj)<-str_replace(sampleNames(rds.obj),"[)]","_")
rds.pdata <- pData(rds.obj[[1]])

#set parameters for fixed bin size
bin_size <- bin*1000
bins<-getBinAnnotations(binSize = bin)
nbins_ref_genome <- sum(fData(rds.obj[[1]])$use)
nbins<-nrow(bins)

#define helper functions
getSegTable<-function(x) #returns a table containing copy number segments
{
    dat<-x
    sn<-assayDataElement(dat,"segmented")
    fd <- fData(dat)
    fd$use -> use
    fdfiltfull<-fd[use,]
    sn<-sn[use,]
    segTable<-c()
    for(c in unique(fdfiltfull$chromosome))
    {
        snfilt<-sn[fdfiltfull$chromosome==c]
        fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
        sn.rle<-rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(1:length(sn.rle$lengths), function(s) {
            from <- fdfilt$start[starts[s]]
            to <- fdfilt$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}

getPloidy<-function(abs_profiles) #returns the ploidy of a sample from segTab or QDNAseq object
{
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segLen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
        ploidy<-sum((segLen/sum(segLen))*as.numeric(segTab$segVal))
        out<-c(out,ploidy)
    }
    data.frame(out,stringsAsFactors = F)
}

getSampNames<-function(abs_profiles) # convenience function for getting sample names from QDNAseq or segTab list
{
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
        samps<-colnames(abs_profiles)
    }
    else
    {
        samps<-names(abs_profiles)
    }
    samps
}

depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

cntodepth<-function(cn,purity,seqdepth) #converts copy number to read depth given purity and single copy depth
{
    seqdepth*((1-purity)*2+purity*cn)
}

get_tp53_bin = function(bin_annotations, qdnaseq_object) {
	# A function to determine the appropriate bin for TP53 in the human reference genome given a variable bin size
	# input:
	#	bin_annotations: genomic bin annotations for a given genome and bin size
	#	qdnaseq_object: defines which to use, which is needed for bin indexing	
	# output: The appropriate bin number for TP53 for the given bin annotations and qdnaseq object

	# load in our table of bin co-ordinates. 'merged' refers to the fact that co-ordinate data and
	# metadata is merged into a single column/variable and needs to be split
	# Only used bins are processed in this step 
	bin_coordinates = tibble::tibble(merged=rownames(bin_annotations[fData(qdnaseq_object[[1]])$use]))

	# reformat bin co-ordinate data so 'merged' data column is split into meaningful constitutents
	bin_coordinates = bin_coordinates %>% dplyr::mutate(
		chromosome = stringr::str_split_fixed(merged, pattern=':', n=2)[,1],
		locus = stringr::str_split_fixed(merged, pattern=':', n=2)[,2]	
	) %>% dplyr::mutate(
		start = stringr::str_split_fixed(locus, pattern='-', n=2)[,1],
		stop = stringr::str_split_fixed(locus, pattern='-', n=2)[,2],
		row_number = dplyr::row_number() # index each bin
	)
	# only select relevant columns
	bin_coordinates = bin_coordinates %>% dplyr::select(row_number,chromosome,start,stop)

	# derive the bin size from the precomputed bin annotations
	bin_size = (bin_coordinates$stop[1] %>% as.integer - bin_coordinates$start[1] %>% as.integer) + 1	

	# Data regarding the TP53 genomic co-ordinates hard-coded from Ensembl/Gencode
	tp53 = tibble::tibble(chromosome=17, build='GRCh37', gencode_release_version=19, start=7565097, stop=7590856)
	
	# perform integer division of TP53 start with the bin size in order to determine the correct bin
	# multiply the result of the integer division by the bin size to map the chosen bin back into genomic co-ordinates
	tp53_segment_lower = ( tp53$start %/% bin_size ) * bin_size
	tp53_segment_upper = tp53_segment_lower + bin_size 
	tp53_segment_lower = tp53_segment_lower + 1 # this offset is needed by our bin reference table

	# filter for TP53 bin
	bin_coordinates = bin_coordinates %>% dplyr::filter(chromosome == tp53$chromosome[1]) %>%
				dplyr::filter(start==tp53_segment_lower) %>%
				dplyr::filter(stop==tp53_segment_upper)
 
	# return bin number
	return (bin_coordinates %>% .$row_number)
} 

# determine the correct tp53 bin index given the chosen bin size
tp53_bin_index = get_tp53_bin(bins,rds.obj)

#estimate absolute copy number fits for all samples in parallel
ploidies<-seq(1.6,8,0.1)
purities<-seq(0.05,1,0.01)
clonality<-c()
#ind<-which(colnames(rds.obj)==sample)
relcn<-rds.obj[[1]]
# added by PS
to_use <- fData(relcn)$use #
relcn <- relcn[to_use,] #
# added by PS
copynumber<-assayDataElement(relcn,"copynumber")
rel_ploidy<-mean(copynumber,na.rm=T)
num_reads<-sum(copynumber,na.rm=T)
sample <- sampleNames(relcn)
print(sample)

res<-foreach(i=1:length(ploidies),.combine=rbind) %do%
{
        ploidy<-ploidies[i]
        print(1)
        rowres<-foreach(j=1:length(purities),.combine=rbind)%do%
        {
            purity<-purities[j]
            downsample_depth<-(((2*(1-purity)+purity*ploidy)/(ploidy*purity))/purity)*15*(2*(1-purity)+purity*ploidy)*nbins_ref_genome*(1/0.91)
            cellploidy<-ploidy*purity+2*(1-purity)
            seqdepth<-rel_ploidy/cellploidy
            cn<-assayDataElement(relcn,"copynumber")
            seg<-assayDataElement(relcn,"segmented")
            cn<-as.numeric(cn[!is.na(cn),])
            seg<-as.numeric(seg[!is.na(seg),])
            integer_cn<-round(depthtocn(seg,purity,seqdepth))
            abs_cn<-depthtocn(seg,purity,seqdepth)
            diffs<-abs(abs_cn-integer_cn)
            TP53cn<-round(depthtocn(seg[tp53_bin_index],purity,seqdepth),1) # to 1 decimal place
            expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))
            clonality<-mean(diffs)
            c(ploidy,purity,clonality,downsample_depth,downsample_depth < rds.pdata$total.reads[row.names(rds.pdata)==sample],TP53cn,expected_TP53_AF)
        }
        rowres
}

colnames(res)<-c("ploidy","purity","clonality","downsample_depth","powered","TP53cn","expected_TP53_AF")
rownames(res)<-1:nrow(res)
res<-data.frame(res,stringsAsFactors = F)
res<-data.frame(apply(res,2,as.numeric,stringsAsFactors=F))
res<-res[order(res$clonality,decreasing=FALSE),]

#output plot of clonality error landscape
pdf(snakemake@output[["pdf"]])
print(ggplot(res,aes(x=ploidy,y=purity,fill=clonality))+geom_tile()+
               scale_fill_gradient(low = "blue", high = "white",trans="log10")+
               theme_bw())
dev.off()

write.table(res,snakemake@output[["csv"]],sep="\t",quote=F,row.names=FALSE)
