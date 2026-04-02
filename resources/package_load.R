# Check if packages are installed
listOfAllPackages = c("Biobase","tidyverse","QDNAseq",
	"QDNAseqmod","QDNAseq.hg19","QDNAseq.hg38",
	"matrixStats","parallel","plyr","ranger","parsnip")

for(thisPackage in listOfAllPackages) {
  
  if(thisPackage %in% rownames(installed.packages()) == FALSE) {
	cat(paste("[install_extra_env][check_packages] ERROR - Package", thisPackage, "needs installing.\n"))
  } else {
	cat(paste("[install_extra_env][check_packages] Package", thisPackage, "is installed.\n"))
  }
}
#a <- lapply(listOfAllPackages,FUN = function(x){suppressPackageStartupMessages(require(x, character.only = TRUE))})

if(!is.null(warnings())){
	print(warnings())
}
