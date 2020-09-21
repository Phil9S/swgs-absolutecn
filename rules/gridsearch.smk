rule gridsearch_fitting:
    input:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/relative_cn_rds/{{project}}_{{sample}}_{{bin}}kb_relSmoothedCN.rds")
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_{sample}_clonality.csv"
    params:
        bin=config["bins"],
        outdir=OUT_DIR,
        project=config["project_name"]
    threads: 1
    shell:
        "Rscript scripts/ploidy_purity_search_standard_error.R {threads} {input} {params.bin} {params.outdir} {params.project} {wildcards.sample}"
