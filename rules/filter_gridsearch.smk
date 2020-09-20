rule gridsearch_filter:
    input:
        cl=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/clonality_results/{project}_clonality.csv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{bin}kb_relSmoothedCN.rds"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    params:
        bin="{bin}",
        meta=config["samplesheet"],
        project="{project}",
        outdir=OUT_DIR
    threads: 1
    shell:  
        "Rscript scripts/gridsearch_results_filtering.R {input.cl} {input.rds} {params.meta} {params.bin} {params.outdir} {params.project}"
        
