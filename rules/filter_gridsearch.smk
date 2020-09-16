rule gridsearch_filter:
    input:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_clonality.csv"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    params:
        bin=config["bins"]
    threads: config["resources"]["threads"]
    shell:
        "touch {output}"
        
