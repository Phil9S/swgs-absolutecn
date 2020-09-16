from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os.path

min_version("5.10.0")
# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

#Predefine output folders (TEMP)
OUT_DIR=config["out_dir"]

#Load sample sheet and set index
samplesheet = pd.read_table(config["samplesheet"]).set_index(["SAMPLE_ID"], drop=False)
validate(samplesheet, schema="../schemas/samples.schema.yaml")
SAMPLES = list(samplesheet['SAMPLE_ID'].unique())

