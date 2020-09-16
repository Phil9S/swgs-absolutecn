rule gridsearch_fitting:
    input:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_clonality.csv",
    params:
        bin=config["bins"]
    threads: config["resources"]["threads"]
    shell:
        "touch {output}"
