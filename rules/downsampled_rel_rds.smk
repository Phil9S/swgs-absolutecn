rule ds_relRDS:
    input:
        bam=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/downsampled_bams/{{sample}}.bam"),
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    threads: 1
    params:
        outdir=OUT_DIR,
        project="{project}",
        bin="{bin}"
    shell:
        "Rscript scripts/qdnaseq_mod_ds.R {params.bin} {threads} {params.outdir} {params.project} {input.meta} {input.bam} {wildcards.sample}"
