rule relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/full_rd/bams/{sample}.bam",sample=SAMPLES)
    output:
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    params:
        bin=config["bins"]
    shell:
        "ls {input}"
        "touch {output}"
