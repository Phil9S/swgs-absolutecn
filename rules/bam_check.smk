rule bam_check:
    input:
        bam=get_bam,
    output:
        ok=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"),
    log:
        "logs/bam_check_{project}_{bin}kb_{sample}.log",
    container:
        image
    threads: 1
    params:
        filetype=config["filetype"],
    script:
        "../scripts/bamCheck.R"
