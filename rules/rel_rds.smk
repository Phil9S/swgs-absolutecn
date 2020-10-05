rule relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    params:
        bin="{bin}",
        outdir=OUT_DIR,
        project="{project}",
        meta=config["samplesheet"]
    script:
        "../scripts/qdnaseq_mod.R"

