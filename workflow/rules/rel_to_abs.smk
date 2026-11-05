rule rel_to_abs:
    input:
        rds=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"),
        meta=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/{project}_{sample}_fitted_QC.tsv"),
    output:
        tsv=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_abs_fits.tsv"),
        rds=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_absCopyNumber.rds"),
        tab=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/{project}_{sample}_{bin}kb_ds_absCopyNumber_segTable.tsv"),
        plot=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/abs_cn_rds/plots/{sample}.png"),
    log:
        "logs/rel_to_abs_{project}_{bin}kb_{sample}.log",
    container:
        image
    params:
        sample="{sample}",
        genome=config["genome"],
    script:
        "../scripts/qdnaseqRelToAbs.R"
