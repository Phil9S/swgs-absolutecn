rule split_input:
    input:
       tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}kb_combined_QC_predownsample.tsv",
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{sample}_fit_QC_predownsample.tsv",
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        bin="{bin}",
        outdir=OUT_DIR,
        project="{project}"
    threads: 1
    script:
        "../scripts/split_input.R"
