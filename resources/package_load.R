# Check if packages are installed
listOfAllPackages = c("Biobase","dplyr","tidyverse","QDNAseq","doMC","QDNAseq",
	"QDNAseqmod","QDNAseq.hg19","QDNAseq.hg38","foreach","ggplot2","stringr",
	"matrixStats","parallel","plyr")

for(thisPackage in listOfAllPackages) {
  
  if(thisPackage %in% rownames(installed.packages()) == FALSE) {
	cat(paste("[install_env] ERROR - Package", thisPackage, "needs installing.\n"))
  } else {
	cat(paste("[install_env] Package", thisPackage, "is installed.\n"))
  }
}
#a <- lapply(listOfAllPackages,FUN = function(x){suppressPackageStartupMessages(require(x, character.only = TRUE))})

if(!is.null(warnings())){
	print(warnings())
}
