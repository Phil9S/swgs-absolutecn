# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
include: "rules/common.smk"

rule all:
    input:
        expand(OUT_DIR+"relative_cn_rds/{bin}kb/{project}_{bin}kb_relSmoothedCN.rds",bin=config["bins"],project=config["project_name"])
        # Subsequent target rules can be specified below. They should start with all_*.

# pipeline rules
include: "rules/rel_rds.smk"
include: "rules/collapse_rel_rds.smk"
include: "rules/other.smk"
