#!/bin/bash

set -e 

script="install_env"
# VARS
INSTALL_BIN="mamba"
CONDA_VERSION=4.8.2
MICROMAMBA_VERSION=1.3.1
CONDA_VERSION_N=$(sed 's/\.//g' <<< "${CONDA_VERSION}")
MICROMAMBA_VERSION_N=$(sed 's/\.//g' <<< "${MICROMAMBA_VERSION}")
INSTALLED_CONDA_VERSION=$(conda -V | sed 's/conda //' | sed 's/\.//g')
INSTALLED_MICROMAMBA_VERSION=$(micromamba --version | sed 's/conda //' | sed 's/\.//g')

## Default behaviour
if [[ $# -eq 0 ]]; then
	echo -e "[${script}] No arguments given. Specify a environment bin; mamba or conda."
	exit 1
fi

if [ $1 == "mamba" ]; then
    # Check MICROMAMBA available
    if ! [ -x "$(command -v micromamba)" ]; then
        echo -e "[${script}] Error: micromamba has not been installed or is not available on \$PATH"
        exit 1
    fi

    # Check micromamba version (rudamentary)
    if [ "${INSTALLED_MICROMAMBA_VERSION}" -lt "${MICROMAMBA_VERSION_N}" ]; then
        echo -e "[${script}] Error - micromamba is older than the required version"
        echo -e "[${script}] Required: ${MICROMAMBA_VERSION} / Installed: $(micromamba --version)"
        exit 1
    fi

    echo -e "[${script}] Creating env"
    # micromamba install
    micromamba env create -y -f config/conda.yaml 1> /dev/null
    eval "$(micromamba shell hook --shell=bash)"
    micromamba activate swgs-abscn

elif [ $1 == "conda" ]; then
    # Check conda available
    if ! [ -x "$(command -v conda)" ]; then
        echo -e "[${script}] Error: conda has not been installed or is not available on PATH"
        exit 1
    fi

    # Check conda version (rudamentary)
    if [ "${INSTALLED_CONDA_VERSION}" -lt "${CONDA_VERSION_N}" ]; then
	echo -e "[${script}] Error - conda/miniconda is older than the required version"
	echo -e "[${script}] Required: ${CONDA_VERSION} / Installed: $(conda -V)"
	exit 1
    fi

    # Check provided conda directory
    if [ "$#" -lt 2 ]; then
        echo -e "[${script}] Error - conda/miniconda directory missing"
        echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
        exit 1
    fi

    # conda install
    # Set conda directory
    if ! [ -d "$2" ]; then
        echo -e "[${script}] Error - conda/miniconda directory not correct"
        echo -e "[${script}] Usage example './install_env.sh conda /home/user/miniconda3/'"
        exit 1
    else 
        CONDA_DIR=$2
    fi
    conda env create -f config/conda.yaml 1> /dev/null
    DIR=${CONDA_DIR}etc/profile.d/conda.sh
    if [ -f "${DIR}" ]; then
        echo -e "[${script}] Initialising conda env"
        source ${CONDA_DIR}etc/profile.d/conda.sh
    else
        echo -e "[${script}] Error: Unable to find conda intialisation script"
        echo -e "[${script}] Error: Make sure to provide the miniconda directory - e.g. '/home/user/miniconda3/'"
        echo -e "[${script}] Usage example './install_env.sh /home/user/miniconda3/'"
        exit 1
    fi
    conda activate swgs-abscn

else
    echo -e "[${script}] Error - neither mamba or conda specified"
    exit 1
fi

echo -e "[${script}] Installing QDNAseq.hg38"
Rscript -e 'remotes::install_github(repo = "asntech/QDNAseq.hg38",quiet=FALSE,upgrade=FALSE,force=TRUE)'
echo -e "[${script}] Installing modified QDNAseq package"
Rscript -e 'remotes::install_github(repo = "markowetzlab/QDNAseqmod",quiet=FALSE,upgrade=FALSE,force=TRUE)'


echo -e "[${script}] Testing package installation"
Rscript resources/package_load.R
echo -e "[${script}] env ready and all packages installed!"
if [ $1 == "mamba" ]; then
    echo -e "[${script}] activate with 'micromamba activate swgs-abscn'"
else
    echo -e "[${script}] activate with 'conda activate swgs-abscn'"
fi

# END
