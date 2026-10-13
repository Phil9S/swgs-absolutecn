from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os
import subprocess
from snakemake.exceptions import WorkflowError

min_version("8.0.0")

##### load config and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

# Predefine output folders
OUT_DIR = config["out_dir"]
OUT_DIR = os.path.join(OUT_DIR, "")

# Load sample sheet and set index
samplesheet = pd.read_table(
    config["samplesheet"],
    dtype={"PATIENT_ID": str, "SAMPLE_ID": str, "TP53freq": float},
).set_index(["SAMPLE_ID"], drop=False)
validate(samplesheet, schema="../schemas/samples.schema.yaml")

if any(samplesheet.SAMPLE_ID.duplicated()):
    sys.exit("Sample sheet contains duplicated sample ids")

# set container uri
image = "docker://" + config["workflowimage"]
cinsigimage = config["cinsigimage"]

##### CHECK MAX > MIN #####
PLMIN = config["ploidy_min"]
PLMAX = config["ploidy_max"]
PUMIN = config["purity_min"]
PUMAX = config["purity_max"]

# Set computeSigs and check
SIGS = config["computeSigs"]

from snakemake.exceptions import WorkflowError

def _active_venv():
    return (
        os.environ.get("PIXI_ENVIRONMENT_NAME")
        or os.environ.get("CONDA_PREFIX")
        or os.environ.get("VIRTUAL_ENV")
    )

if SIGS == "TRUE":
    venv = _active_venv()
    cinsig_sif = config.get("cinsigimage", "")

    if venv:
        # Local/conda/pixi mode -- verify the package is actually importable
        try:
            subprocess.run(
                ["Rscript", "-e",
                 'if (!requireNamespace("CINSignatureQuantification", quietly=TRUE)) quit(status=1)'],
                check=True, capture_output=True
            )
        except subprocess.CalledProcessError:
            raise WorkflowError(
                f"computeSigs=True and environment '{venv}' is active, but "
                "CINSignatureQuantification is not installed in it.\n"
                "See README section 'Optional - CINSignatureQuantification' "
                "for install instructions\n"
		"Note: GAP Available Source License v1.0, non-commercial "
                "academic use only."
            )
        except FileNotFoundError:
            raise WorkflowError(
                f"computeSigs=True and environment '{venv}' is active, but "
                "Rscript was not found on PATH within it."
            )

    elif cinsig_sif:
        # No active env, but a container image is configured -- verify it exists
        if not os.path.exists(cinsig_sif):
            raise WorkflowError(
                f"computeSigs=True and 'cinsigimager' is set to "
                f"'{cinsig_sif}', but that file does not exist.\n"
                "Build it with:\n"
                f"  singularity build {cinsig_sif} docker/Dockerfile.cinsig-asl\n"
                "Note: GAP Available Source License v1.0, non-commercial "
                "academic use only."
            )

    else:
        # Neither signal present -- can't tell how the user intends to run this
        raise WorkflowError(
            "computeSigs=True but no way to run CINSignatureQuantification was "
            "detected.\n"
            "Either:\n"
            "  (a) activate the environment containing it before running "
            "snakemake, e.g. 'pixi shell', or\n"
            "  (b) set 'cinsigimage' in config to the cinsignaturequantification container\n"
            "See README section 'Optional - CINSignatureQuantification'."
        )

if PLMIN > PLMAX:
    sys.exit("Config error - Minimum ploidy exceeds or is equal to maximum ploidy")

if PUMIN > PUMAX:
    sys.exit("Config error - Minimum purity exceeds or is equal to maximum purity")


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
        f = pd.read_table(fl).set_index(["SAMPLE_ID"], drop=False)
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


