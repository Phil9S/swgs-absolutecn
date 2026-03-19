rule combine:
    input:
       expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{{project}}_{sample}_fit_QC_predownsample.tsv",sample=SAMPLES)
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    threads: 1
    script:
        "../scripts/combine.R"
