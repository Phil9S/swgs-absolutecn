rule rel_to_abs:
    input:
        rds=expand(
            OUT_DIR
            + "sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{{sample}}/relative_cn_rds/{{project}}_{{sample}}_{{bin}}kb_relSmoothedCN.rds"
        ),
        meta=expand(
            OUT_DIR
            + "sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{{sample}}/{{project}}_{{sample}}_fitted_QC.tsv"
        ),
    output:
        tsv=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_absCopyNumber.rds",
        tab=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_absCopyNumber_segTable.tsv",
        plot=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/plots/{sample}.png",
    log:
        "logs/rel_to_rds_{project}_{bin}kb_{sample}.log",
    container:
        image
    params:
        sample="{sample}",
        project="{project}",
        bin="{bin}",
        genome=config["genome"],
    script:
        "../scripts/qdnaseqRelToAbs.R"
