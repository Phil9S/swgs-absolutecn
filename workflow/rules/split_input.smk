rule split_input:
    input:
        tsv=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}kb_combined_QC_predownsample.tsv",
    output:
        tsv=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/{project}_{sample}_fitted_QC.tsv",
    log:
        "logs/split_input_{project}_{bin}kb_{sample}.log",
    container:
        image
    threads: 1
    params:
        sample="{sample}",
    script:
        "../scripts/splitInput.R"
