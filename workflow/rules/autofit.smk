rule autofit:
    input:
        tsv=expand(
            OUT_DIR
            + "sWGS_fitting/{{project}}_{{bin}}kb/absolute_PRE_down_sampling/{{sample}}/{{project}}_{{sample}}_fit_QC_predownsample.tsv"
        ),
    output:
        tsv=OUT_DIR
        + "sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/{project}_{sample}_fitted_QC.tsv",
    log:
        "logs/autofit_{project}_{bin}kb_{sample}.log",
    container:
        image
    threads: 1
    params:
        flagThreshold=config["flagThreshold"],
        fitMethod=config["fitMethod"],
        errorMetric=config["errorMetric"],
    script:
        "../scripts/autofit.R"
