rule collapse_relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/full_rd/relative_cn_rds/{{bin}}kb/{sample}_{{bin}}kb.rds",sample=SAMPLES)
    output:
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/relative_cn_rds/{bin}kb/{project}_{bin}kb_relSmoothedCN.rds"
    params:
        bin=config["bins"]
    threads: config["resources"]["threads"]
    shell:
        "ls {input}"
        "touch {output.rds}"

