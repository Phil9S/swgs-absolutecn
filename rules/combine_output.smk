rule combine_output:
    input:
       tsv=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{sample}/{{project}}_{sample}_fit_QC_predownsample.tsv",sample=SAMPLES),
       rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{sample}/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES)
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}kb_combined_QC_predownsample.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    log:
        "logs/combine_output_{project}_{bin}kb.log"
    container:
        image
    threads: 1
    script:
        "../scripts/combine_results.R"
