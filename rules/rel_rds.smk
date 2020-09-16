rule relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/bams/{sample}.bam",sample=SAMPLES)
    output:
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    params:
        bin=config["bins"]
    shell:
        "touch {output}"
