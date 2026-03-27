rule downsample:
    input:
        bam=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam",
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/{project}_{sample}_fitted_QC.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    output:
        ds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/downsampled_bams/{sample}.bam"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        prplpu=prplpu,
        sample="{sample}",
        reference=config["reference"],
        filetype=config["filetype"]
    script:
        "../scripts/downsampleBams.R"
    
