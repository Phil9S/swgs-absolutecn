# Check if packages are installed
listOfAllPackages = c("Biobase","tidyverse","QDNAseq",
	"QDNAseqmod","QDNAseq.hg19","QDNAseq.hg38","rswgsabsolutecn",
	"matrixStats","parallel","plyr","ranger","parsnip")

if (!requireNamespace("CINSignatureQuantification", quietly = TRUE)) {
  message(
    "\nCINSignatureQuantification is required when computeSigs = TRUE.\n\n",
    "This package is licensed under the GAP Available Source License v1.0 ",
    "(non-commercial academic use only,\n see https://github.com/markowetzlab/CINSignatureQuantification/blob/main/LICENSE)",
    "and is NOT bundled with this pipeline for that reason.\n\n",
    "To install it yourself, run:\n",
    "Rscript -e 'remotes::install_github(\"markowetzlab/CINSignatureQuantification\")'\n\n",
    "By installing it you agree to its license terms, including the non-commercial restriction."
    )
}

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

if (!requireNamespace("CINSignatureQuantification", quietly = TRUE)) {
  message(
    "\nCINSignatureQuantification is required when computeSigs = TRUE.\n\n",
    "This package is licensed under the GAP Available Source License v1.0 ",
    "(non-commercial academic use only,\n see https://github.com/markowetzlab/CINSignatureQuantification/blob/main/LICENSE)",
    "and is NOT bundled with this pipeline for that reason.\n\n",
    "To install it run:\n",
    "Rscript -e 'remotes::install_github(\"markowetzlab/CINSignatureQuantification\")'\n\n",
    "By installing it you agree to its license terms, including the non-commercial restriction."
    )
}

