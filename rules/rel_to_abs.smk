rule rel_to_abs:
    input:
        rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES),
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_absCopyNumber.rds"
    threads: 1
    params:
        outdir=OUT_DIR,
        project="{project}",
        bin="{bin}"
    shell:
        "Rscript scripts/qdnaseq_rel_to_abs.R {input.meta} {params.outdir} {params.bin} {params.project}"

