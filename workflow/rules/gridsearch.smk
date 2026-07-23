rule gridsearch:
    input:
        rds=expand(
            OUT_DIR
            + "sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{{sample}}/relative_cn_rds/{{project}}_{{sample}}_{{bin}}kb_relSmoothedCN.rds"
        ),
    output:
        tsv=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/clonality_results/{project}_{sample}_clonality.tsv",
        pdf=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/clonality_results/{project}_{sample}_clonality.pdf",
        filt=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/{project}_{sample}_filtered_results.tsv",
        fit=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/{project}_{sample}_fit_QC_predownsample.tsv",
        plot=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/plots/{sample}.png",
    log:
        "logs/gridsearch_fitting_{project}_{bin}kb_{sample}.log",
    container:
        image
    params:
        bin="{bin}",
        project="{project}",
        meta=config["samplesheet"],
        ploidy_min=config["ploidy_min"],
        ploidy_max=config["ploidy_max"],
        purity_min=config["purity_min"],
        purity_max=config["purity_max"],
        homozygous_threshold=config["homozygous_threshold"],
        genome=config["genome"],
        sample="{sample}",
        af_cutoff=config["af_cutoff"],
        filter_underpowered=config["filter_underpowered"],
        filter_homozygous=config["filter_homozygous"],
        homozygous_prop=config["homozygous_prop"],
    script:
        "../scripts/gridsearch.R"
