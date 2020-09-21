#!/bin/bash

echo -e "Creating conda env"
conda env create -f config/conda.yaml
echo -e "Activating conda env"
conda activate swgs-abscn
echo -e "Testing package loading"


