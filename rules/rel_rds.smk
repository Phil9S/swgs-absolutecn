rule relRDS:
    input:
        bams=expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam"),
    output:
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    params:
        bin="{bin}",
        meta=config["samplesheet"]
        pairedEnd=config["pairedEnd"],
        use_seed=config["use_seed"],
        seed_val=config["seed_val"],
        genome=config["genome"],
        sample="{sample}",
        autoSmooth=config["autoSmooth"],
        smoothThreshold=config["smoothThreshold"]
    script:
        "../scripts/qdnaseq_mod.R"

