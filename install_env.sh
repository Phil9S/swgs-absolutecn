#!/bin/bash

set -e 

script="install_env"
# VARS
echo -e "[${script}] Installing additional packages - QDNAseq.hg38"
Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseq.hg38",quiet=TRUE,upgrade=FALSE,force=FALSE)'
echo -e "[${script}] Installing additional packages - QDNAseqmod package"
Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseqmod",quiet=TRUE,upgrade=FALSE,force=FALSE)'
echo -e "[${script}] Installing additional packages - rswgsabsolutecn package"
Rscript -e 'remotes::install_github(repo = "markowetzlab/r-swgs-absolutecn",quiet=TRUE,upgrade=FALSE,force=FALSE)'
echo -e "[${script}] Testing package installation..."
Rscript resources/package_load.R
# END
