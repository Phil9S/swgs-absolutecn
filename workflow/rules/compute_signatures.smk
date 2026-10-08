rule compute_signatures:
    input:
        tsv=expand(
            OUT_DIR
            + "sWGS_fitting/{{project}}_{{bin}}kb/absolute_POST_down_sampling/{{sample}}/abs_cn_rds/{{project}}_{{sample}}_{{bin}}kb_ds_absCopyNumber_segTable.tsv"
        ),
    output:
        tsvDrews=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/cin_signatures/{sample}/{project}_{sample}_{bin}kb_DREWS_CIN_signature_acts.tsv",
        tsvMac=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/cin_signatures/{sample}/{project}_{sample}_{bin}kb_MAC_CIN_signature_acts.tsv",
    log:
        "logs/compute_signatures_{project}_{bin}kb_{sample}.log",
    container:
        ""
    params:
        bin="{bin}",
        genome=config["genome"],
        project="{project}",
        sample="{sample}",
    script:
        "../scripts/computeSignatures.R"

