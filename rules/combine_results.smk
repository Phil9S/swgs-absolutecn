rule combine_results:
    input:
       tsv=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{{project}}_{sample}_{{bin}}kb_ds_abs_fits.tsv",sample=SAMPLES),
       rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{{project}}_{sample}_{{bin}}kb_ds_absCopyNumber.rds",sample=SAMPLES),
       tab=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{{project}}_{sample}_{{bin}}kb_ds_absCopyNumber_segTable.tsv",sample=SAMPLES)
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{project}_{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{project}_{bin}kb_ds_absCopyNumber.rds",
        tab=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{project}_{bin}kb_ds_absCopyNumber_segTable.tsv"
    log:
        "logs/combine_output_{project}_{bin}kb.log"
    container:
        image
    threads: 1
    script:
        "../scripts/combineResults.R"
