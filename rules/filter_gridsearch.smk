rule gridsearch_filter:
    input:
        cl=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/clonality_results/{{project}}_{sample}_clonality.tsv",sample=SAMPLES),
        rds=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/relative_cn_rds/{{project}}_{sample}_{{bin}}kb_relSmoothedCN.rds",sample=SAMPLES)
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        bin="{bin}",
        meta=config["samplesheet"],
        project="{project}",
        outdir=OUT_DIR,
        af_cutoff=config["af_cutoff"],
        filter_underpowered=config["filter_underpowered"],
        filter_homozygous=config["filter_homozygous"],
        homozygous_prop=config["homozygous_prop"]
    threads: THREADS 
    script: 
        "../scripts/gridsearch_results_filtering.R"
        
