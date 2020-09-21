rule relRDS:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    params:
        bin=config["bins"],
        outdir=OUT_DIR,
        project=config["project_name"],
        meta=config["samplesheet"]
    threads: 1
    shell:
        "Rscript scripts/qdnaseq_mod.R {params.bin} {threads} {params.outdir} {params.project} {params.meta} {input} {wildcards.sample}"

