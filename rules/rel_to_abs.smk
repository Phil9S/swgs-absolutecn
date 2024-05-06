rule rel_to_abs:
    input:
        rds=lambda wildcards: expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLE_LISTS[wildcards.bin]),
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_absCopyNumber.rds"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        outdir=OUT_DIR,
        project="{project}",
        bin="{bin}"
    threads: THREADS 
    script:
        "../scripts/qdnaseq_rel_to_abs.R"

