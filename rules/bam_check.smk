rule bam_check:
    input:
       bam=FILE_LIST
    output:
        ok=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"
    log:
        "logs/bam_check_{project}_{bin}kb_{sample}.log"
    container:
        image
    threads: 1
    script:
        "../scripts/bam_check.R"
