rule gridsearch_fitting:
    input:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/clonality_results/{project}_clonality.csv",
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/clonality_results/{project}_clonality.pdf"
    params:
        bin=config["bins"]
    threads: config["resources"]["threads"]
    shell:
        "ls {input}"
        "touch {output}"
