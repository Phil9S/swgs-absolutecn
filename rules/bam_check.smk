rule check_bam:
    input:
       bam=FILE_LIST
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bam.ok"
    threads: 1
    script:
        "../scripts/bam_check.R"
