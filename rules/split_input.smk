rule split_input:
    input:
       tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}kb_combined_QC_predownsample.tsv",
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{project}_{sample}_fitted_QC.tsv"
    container:
        image
    params:
	sample="{sample}"
    threads: 1
    script:
        "../scripts/split_input.R"
