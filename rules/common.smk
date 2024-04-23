from snakemake.utils import validate
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
import os.path
import ruamel.yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString, DoubleQuotedScalarString

yaml = ruamel.yaml.YAML()

min_version("5.10.0")

#### Check max threads and job submission type ####
patterns = ['local','slurm','pbs']
strings = sys.argv

def get_job_sub(strings,patterns):
    for p in patterns:
        if any([i for i in strings if p in i]):
            if p == 'local':
                return p
            elif p == 'slurm':
                return p
            elif p == 'pbs':
                return p

def get_threads():
    job_sub_val = get_job_sub(strings,patterns)
    if job_sub_val is None:
        file_name = "profile/local/cluster_config.yaml"
        cluster_config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
        return cluster_config['cpus-per-task']
    else:
        if job_sub_val == 'local':
            file_name = "profile/local/cluster_config.yaml"
            cluster_config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
            return cluster_config['cpus-per-task']
        elif job_sub_val == 'slurm':
            file_name = "profile/slurm/cluster_config.yaml"
            cluster_config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
            return cluster_config['rel_to_abs']['cpus-per-task']
        elif job_sub_val == 'pbs':
            file_name = "profile/pbs/cluster_config.yaml"
            cluster_config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open(file_name))
            resourceSplit = cluster_config['rel_to_abs']['l'].split(',')
            proc = resourceSplit[0].split(':')
            proc_t = proc[1].split('=')
            return(int(proc_t[1]))

THREADS = get_threads()
##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

#Predefine output folders
OUT_DIR=config["out_dir"]
OUT_DIR=os.path.join(OUT_DIR,"")

#Load sample sheet and set index
samplesheet = pd.read_table(config["samplesheet"],dtype={'PATIENT_ID': str,'SAMPLE_ID':str,'TP53freq':float}).set_index(["SAMPLE_ID"], drop=False)
validate(samplesheet, schema="../schemas/samples.schema.yaml")

#### Check bin values ####

BIN_VALS = config["bins"]
BIN_DEF = [1,5,15,30,50,100,500,1000]

if not set(BIN_VALS).issubset(BIN_DEF):
    sys.exit("Config error - Some specified bin values are not available")

##### CHECK MAX > MIN #####
PLMIN=config["ploidy_min"]
PLMAX=config["ploidy_max"]
PUMIN=config["purity_min"]
PUMAX=config["purity_max"]

if PLMIN > PLMAX:
    sys.exit("Config error - Minimum ploidy exceeds or is equal to maximum ploidy")

if PUMIN > PUMAX:
    sys.exit("Config error - Minimum purity exceeds or is equal to maximum purity")

