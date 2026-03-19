rule combine_output:
    input:
       tsv=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{{project}}_{sample}_fit_QC_predownsample.tsv",sample=SAMPLES)
       #rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES)
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}_combined_QC_predownsample.tsv"
        #rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        bin="{bin}",
        outdir=OUT_DIR,
        project="{project}",
        sample="{sample}"
    threads: 1
    script:
        "../scripts/combine_output.R"
