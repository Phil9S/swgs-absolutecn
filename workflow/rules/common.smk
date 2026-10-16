from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os
import subprocess
import sys
from snakemake.exceptions import WorkflowError

min_version("8.0.0")

##### load config and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

# Predefine output folders
OUT_DIR = config["out_dir"]
OUT_DIR = os.path.join(OUT_DIR, "")

# Check CRAM reference
if config["filetype"] == "CRAM" and not config["reference"]:
    sys.exit("CRAM input requires 'reference'")

# Load sample sheet and set index
samplesheet = pd.read_table(
    config["samplesheet"],
    dtype={"PATIENT_ID": str, "SAMPLE_ID": str, "TP53freq": float},
).set_index(["SAMPLE_ID"], drop=False)
validate(samplesheet, schema="../schemas/samples.schema.yaml")

if any(samplesheet.SAMPLE_ID.duplicated()):
    sys.exit("Sample sheet contains duplicated sample ids")

# precPloidy/precPurity columns must exist but may contain NA
# snakemake validate() drops NA values so schema 'required' rejects NA rows
required_cols = {"precPloidy", "precPurity"}
missing_cols = required_cols - set(samplesheet.columns)
if missing_cols:
    sys.exit(
        "Sample sheet is missing required column(s): "
        + ", ".join(sorted(missing_cols))
        + " - add the column(s) with NA where values are not precomputed"
    )

# Accept docker paths with or without the docker:// prefix
image = "docker://" + config["workflowimage"]
cinsigimage = config["cinsigimage"]
if cinsigimage and "://" not in cinsigimage:
    cinsigimage = "docker://" + cinsigimage

##### CHECK MAX > MIN #####
PLMIN = config["ploidy_min"]
PLMAX = config["ploidy_max"]
PUMIN = config["purity_min"]
PUMAX = config["purity_max"]

# Set computeSigs and check
SIGS = config["computeSigs"]

# For computeSigs, basic check for CINSignatureQuantification availability
# Env/container run mode is taken from config where setting 'cinsigimage' means container mode
if SIGS == "TRUE" and not cinsigimage:
    try:
        subprocess.run(
            ["Rscript", "-e",
             'if (!requireNamespace("CINSignatureQuantification", quietly=TRUE)) quit(status=1)'],
            check=True, capture_output=True
        )
    except subprocess.CalledProcessError:
        raise WorkflowError(
            "computeSigs=TRUE but CINSignatureQuantification is not installed "
            "in the active environment.\n"
            "Either install it, or set 'cinsigimage' in config to the "
            "CINSignatureQuantification docker image and run with "
            "containers (--sdm apptainer).\n"
            "See README section 'Optional - CINSignatureQuantification'.\n"
            "Note: GAP Available Source License v1.0, non-commercial "
            "academic use only."
        )
    except FileNotFoundError:
        raise WorkflowError(
            "computeSigs=TRUE but Rscript was not found on PATH in the "
            "active environment."
        )

if PLMIN > PLMAX:
    sys.exit("Config error - Minimum ploidy exceeds maximum ploidy")

if PUMIN > PUMAX:
    sys.exit("Config error - Minimum purity exceeds maximum purity")


def get_bam(wc):
    try:
        val = samplesheet.at[wc.sample, "file"]
    except KeyError:
        raise ValueError(
            f"Sample '{wc.sample}' not found in sample sheet ({config['samplesheet']})"
        )
    return str(val)

def get_ds_bams(P, b):
    d = []
    fl = (
        OUT_DIR
        + "sWGS_fitting/"
        + P
        + "_"
        + str(b)
        + "kb/absolute_PRE_down_sampling/"
        + P
        + "_"
        + str(b)
        + "kb_combined_QC_predownsample.tsv"
    )
    if os.path.isfile(fl):
        f = pd.read_table(fl, dtype={"SAMPLE_ID": str, "PATIENT_ID": str}).set_index(["SAMPLE_ID"], drop=False)
        if not f.empty:
            if f["use"].dtype == bool and all(~f.use.isnull()):
                ftb = f[(f["use"] == True)]
                if all(~ftb.SAMPLE_ID.duplicated()):
                    d = list(ftb["SAMPLE_ID"])
                else:
                    sys.exit("Error - QC file contains duplicated sample ids marked 'TRUE'")
            else:
                sys.exit("Error - QC file 'use' column contains non-boolean values")
        else:
            sys.exit("Error - QC file is empty. Filtering criteria may need adjusting")
    else:
        sys.exit("Error - QC file from stage 1 is missing")
    return d


