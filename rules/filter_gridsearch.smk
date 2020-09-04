rule gridsearch_filter:
    input:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/clonality_results/{project}_clonality.csv"
    output:
        filtered=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/{project}_filtered_results.tsv",
        QC=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/full_rd/{project}_fit_QC_predownsample.tsv",
    params:
        bin=config["bins"]
    threads: config["resources"]["threads"]
    shell:
        "ls {input}"
        "ls {output.filtered}"
        
