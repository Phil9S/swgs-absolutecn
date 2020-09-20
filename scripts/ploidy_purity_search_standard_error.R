#this script outputs a ranked list of absolute copy number fits using the 
#traditional clonality error function. it is designed to be used in conjuction 
#with TP53 allele frequency to determine precise purity and ploidy fit

#grab commandline arguments
args = commandArgs(trailingOnly=TRUE)

cores <- as.numeric(args[1])
rds.filename <- args[2]
bin <- as.numeric(args[3])
out_dir <- args[4]
project <- args[5]

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
colnames(rds.obj)<-str_replace(colnames(rds.obj),"[(]","_")
colnames(rds.obj)<-str_replace(colnames(rds.obj),"[)]","_")
rds.pdata <- pData(rds.obj)

#set parameters for fixed bin size
bin_size <- bin*1000
bins<-getBinAnnotations(binSize = bin)
nbins_ref_genome <- sum(fData(rds.obj)$use)
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

#estimate absolute copy number fits for all samples in parallel
foreach(sample=row.names(rds.pdata))%dopar%
{
    ploidies<-seq(1.6,8,0.1)
    purities<-seq(0.05,1,0.01)
    clonality<-c()
    ind<-which(colnames(rds.obj)==sample)
    relcn<-rds.obj[,ind]
    # added by PS
    to_use <- fData(relcn)$use #
    relcn <- relcn[to_use,] #
    # added by PS
    copynumber<-assayDataElement(relcn,"copynumber")
    rel_ploidy<-mean(copynumber,na.rm=T)
    num_reads<-sum(copynumber,na.rm=T)
    #print(sample)
    #print(num_reads)
    
    res<-foreach(i=1:length(ploidies),.combine=rbind)%do%
    {
        ploidy<-ploidies[i]
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
            TP53cn<-round(depthtocn(seg[73504],purity,seqdepth),1) # to 1 decimal place
            expected_TP53_AF<-TP53cn*purity/(TP53cn*purity+2*(1-purity))
            clonality<-mean(diffs)
            c(ploidy,purity,clonality,downsample_depth,downsample_depth<rds.pdata$total.reads[row.names(rds.pdata)==sample],TP53cn,expected_TP53_AF)
        }
        rowres
    }
    
    colnames(res)<-c("ploidy","purity","clonality","downsample_depth","powered","TP53cn","expected_TP53_AF")
    rownames(res)<-1:nrow(res)
    res<-data.frame(res,stringsAsFactors = F)
    res<-data.frame(apply(res,2,as.numeric,stringsAsFactors=F))
    res<-res[order(res$clonality,decreasing=FALSE),]

    #output plot of clonality error landscape
    pdf(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",sample,"_clonality.pdf"))
    print(ggplot(res,aes(x=ploidy,y=purity,fill=clonality))+geom_tile()+
              scale_fill_gradient(low = "blue", high = "white",trans="log10")+
              theme_bw())
    dev.off()
    
    #write table of clonality scores
    print(paste0("Writing clonality table for sample: ",sample))
    write.table(res,paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",sample,"_clonality.csv"),sep="\t",quote=F,row.names=FALSE)
}

filelist <- list.files(pattern="*clonality.csv",path=paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/"))
clonality <- do.call(rbind,
			lapply(filelist,FUN = function(x){
				n <- gsub(pattern="_clonality.csv",rep="",x=x)
				tab <- read.table(paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",x),sep="\t",skip=1)
				tab <- cbind(rep(n,times=nrow(tab)),tab)
				return(tab)
			}))
colnames(clonality) <- c("SAMPLE_ID","ploidy","purity","clonality","downsample_depth","powered","TP53cn","expected_TP53_AF")
write.table(clonality,paste0(out_dir,"sWGS_fitting/",project,"_",bin,"kb/absolute_PRE_down_sampling/clonality_results/",project,"_clonality.csv"),sep="\t",quote=F,row.names=FALSE)
