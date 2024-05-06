rule check_bam:
    input:
       bam=FILE_LIST
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bam.ok"
    singularity:
        image_base_url+"swgs-absolutecn:latest"
    threads: 1
    script:
        "../scripts/bam_check.R"
