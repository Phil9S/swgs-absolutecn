rule downsample:
    input:
        bam=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam",
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{sample}/{project}_{sample}_fitted_QC.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    output:
        bam=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/downsampled_bams/{sample}.bam"
    log:
        "logs/downsample_{project}_{bin}kb_{sample}.log"
    container:
        image
    params:
        prplpu=prplpu,
        sample="{sample}",
        reference=config["reference"],
        filetype=config["filetype"]
    script:
        "../scripts/downsampleBams.R"
    
