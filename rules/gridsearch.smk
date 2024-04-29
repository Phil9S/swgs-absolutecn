rule gridsearch_fitting:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/relative_cn_rds/{{project}}_{{sample}}_{{bin}}kb_relSmoothedCN.rds")
    output:
        tsv=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_{sample}_clonality.tsv",
        pdf=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_{sample}_clonality.pdf"
    params:
        bin="{bin}",
        outdir=OUT_DIR,
        project="{project}",
        meta=config["samplesheet"],
        ploidy_min=config["ploidy_min"],
        ploidy_max=config["ploidy_max"],
        purity_min=config["purity_min"],
        purity_max=config["purity_max"]
    script:
        "../scripts/ploidy_purity_search_standard_error.R"
