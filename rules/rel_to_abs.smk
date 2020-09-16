rule rel_to_abs:
    input:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_abs_fits.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/{project}_{bin}kb_ds_absCopyNumber.rds"
    shell:
        "touch {output.tsv} & touch {output.rds}"

