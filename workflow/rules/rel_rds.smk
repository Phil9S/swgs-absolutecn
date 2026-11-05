rule relRDS:
    input:
        bam=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam"),
        # CRAM-converted BAMs are temp - keep the index until QDNAseq has read the BAM
        bai=(Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam.bai")
          if config["filetype"] == "CRAM" else []),
    output:
        rds=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{sample}/relative_cn_rds/{project}_{sample}_{bin}kb_relSmoothedCN.rds"),
    log:
        "logs/relRDS_{project}_{bin}kb_{sample}.log",
    container:
        image
    params:
        bin="{bin}",
        meta=config["samplesheet"],
        pairedEnd=config["pairedEnd"],
        use_seed=config["use_seed"],
        seed_val=config["seed_val"],
        genome=config["genome"],
        sample="{sample}",
        autoSmooth=config["autoSmooth"],
        smoothThreshold=config["smoothThreshold"],
    script:
        "../scripts/qdnaseqMod.R"
