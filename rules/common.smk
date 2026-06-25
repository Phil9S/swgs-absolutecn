from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os.path

min_version("8.0.0")

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

#Predefine output folders
OUT_DIR=config["out_dir"]
OUT_DIR=os.path.join(OUT_DIR,"")

#Load sample sheet and set index
samplesheet = pd.read_table(config["samplesheet"],dtype={'PATIENT_ID': str,'SAMPLE_ID':str,'TP53freq':float}).set_index(["SAMPLE_ID"], drop=False)
validate(samplesheet, schema="../schemas/samples.schema.yaml")

if any(samplesheet.SAMPLE_ID.duplicated()):
    sys.exit("Sample sheet contains duplicated sample ids")
    
# set container uri
image = config["image"]

##### CHECK MAX > MIN #####
PLMIN=config["ploidy_min"]
PLMAX=config["ploidy_max"]
PUMIN=config["purity_min"]
PUMAX=config["purity_max"]

if PLMIN > PLMAX:
    sys.exit("Config error - Minimum ploidy exceeds or is equal to maximum ploidy")

if PUMIN > PUMAX:
    sys.exit("Config error - Minimum purity exceeds or is equal to maximum purity")

def get_bam(wc):
    try:
        val = samplesheet.at[wc.sample, "file"]
    except KeyError:
        raise ValueError(f"Sample '{wc.sample}' not found in sample sheet ({config['samplesheet']})")
    return str(val)
