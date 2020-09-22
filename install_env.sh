#!/bin/bash

script=$0

if ! [ -x "$(command -v conda)" ]; then
	echo -e "[${script}] Error: conda has not been installed"
	exit 1
fi

if [ "$#" -lt 1 ]; then
	echo -e "[${script}] Error - conda/miniconda directory missing"
	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
	exit 1
fi

CONDA_DIR=$1

echo -e "[${script}] Creating conda env"
conda env create -f config/conda.yaml
DIR=${CONDA_DIR}etc/profile.d/conda.sh
if [ -f "${DIR}" ]; then
	echo -e "[${script}] Initialising conda env"
	source ${CONDA_DIR}etc/profile.d/conda.sh
else
	echo -e "[${script}] Error: Unable to find conda intialisation script"
	echo -e "[${script}] Error: Make sure to provide the miniconda directory - e.g. '/home/user/miniconda3/'"
	echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
	exit
fi

echo -e "[${script}] Activating conda env"
conda activate swgs-abscn
echo -e "[${script}] Adding modified QDNAseq package"
R_LIB_PATH=$(Rscript resources/libpath.R)
cp -r resources/packages/QDNAseqmod/ ${R_LIB_PATH} 
echo -e "[${script}] Testing package installation"
Rscript resources/package_load.R
echo -e "[${script}] conda env ready and all packages installed!"
