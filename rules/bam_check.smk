rule bam_check:
    input:
       bam=get_bam
    output:
        ok=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"
    log:
        "logs/bam_check_{project}_{bin}kb_{sample}.log"
    container:
        image
    params:
        filetype=config["filetype"]
    threads: 1
    script:
        "../scripts/bam_check.R"
