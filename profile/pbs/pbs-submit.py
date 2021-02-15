#!/usr/bin/env python3
"""
Snakemake PBS submit script.
"""
import warnings  # use warnings.warn() rather than print() to output info in this script

from snakemake.utils import read_job_properties

import pbs_utils

CLUSTER_CONFIG = "cluster_config.yaml"

# parse job
jobscript = pbs_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

pbs_options = {}
cluster_config = pbs_utils.load_cluster_config(CLUSTER_CONFIG)

# 2) cluster_config defaults
pbs_options.update(cluster_config["__default__"])

# 4) cluster_config for particular rule
pbs_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 7) Format pattern in snakemake style
pbs_options = pbs_utils.format_values(pbs_options, job_properties)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    pbs_utils.ensure_dirs_exist(pbs_options[o]) if o in pbs_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(pbs_utils.submit_job(jobscript, **pbs_options))
