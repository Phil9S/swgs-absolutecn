#!/bin/bash

set -e 

script="install_env"
# VARS
#CONDA_VERSION=4.8.2
#CONDA_VERSION_N=$(sed 's/\.//g' <<< "${CONDA_VERSION}")
#INSTALLED_CONDA_VERSION=$(conda -V | sed 's/conda //' | sed 's/\.//g')

# Check conda available
#if ! [ -x "$(command -v conda)" ]; then
#	echo -e "[${script}] Error: conda has not been installed or is not available on PATH"
#	exit 1
#fi

# Check conda version (rudamentary)
#if [ "${INSTALLED_CONDA_VERSION}" -lt "${CONDA_VERSION_N}" ]; then
#	echo -e "[${script}] Error - conda/miniconda is older than the required version"
#	echo -e "[${script}] Required: conda ${CONDA_VERSION} / Installed: $(conda -V)"
#	exit	
#fi

# Check provided conda directory
#if [ "$#" -lt 1 ]; then
#	echo -e "[${script}] Error - conda/miniconda directory missing"
#	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
#	exit 1
#fi

echo -e "[${script}] Creating env"

# micromamba install
micromamba env create -y -f config/conda.yaml
eval "$(micromamba shell hook --shell=bash)"
micromamba activate swgs-abscn

# conda install
# Set conda directory
#CONDA_DIR=$1
#conda env create -f config/conda.yaml
#DIR=${CONDA_DIR}etc/profile.d/conda.sh
#if [ -f "${DIR}" ]; then
#	echo -e "[${script}] Initialising conda env"
#	source ${CONDA_DIR}etc/profile.d/conda.sh
#else
#	echo -e "[${script}] Error: Unable to find conda intialisation script"
#	echo -e "[${script}] Error: Make sure to provide the miniconda directory - e.g. '/home/user/miniconda3/'"
#	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
#	exit
#fi
#conda activate swgs-abscn

#R_LIB_PATH=$(Rscript resources/libpath.R)
#cp -r resources/packages/QDNAseqmod/ ${R_LIB_PATH} 

echo -e "[${script}] Installing modified QDNAseq package"
Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseqmod",quiet=TRUE,upgrade=FALSE)'
echo -e "[${script}] Testing package installation"
Rscript resources/package_load.R
echo -e "[${script}] env ready and all packages installed!"
echo -e "[${script}] activate with 'micromamba activate swgs-abscn'"
