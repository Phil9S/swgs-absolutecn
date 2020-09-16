rule ds_relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/downsampled_bams/{sample}.bam",sample=SAMPLES)
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    params:
        bin=config["bins"]
    shell:
        "touch {output}"
