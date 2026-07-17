# Check if packages are installed
listOfAllPackages = c("Biobase","tidyverse","QDNAseq",
	"QDNAseqmod","QDNAseq.hg19","QDNAseq.hg38","rswgsabsolutecn",
	"matrixStats","parallel","plyr","ranger","parsnip")

for(thisPackage in listOfAllPackages) {
  
  if(thisPackage %in% rownames(installed.packages()) == FALSE) {
	cat(paste("[install_env][check_packages] ERROR - Package", thisPackage, "needs installing.\n"))
  } else {
	cat(paste("[install_env][check_packages] Package", thisPackage, "is installed.\n"))
  }
}

if(!is.null(warnings())){
	print(warnings())
}
