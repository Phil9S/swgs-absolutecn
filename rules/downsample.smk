rule downsample:
    input:
        OUT_DIR+"sWGS_fitting/bams/{sample}.bam"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/downsampled_bams/{sample}.bam"
    shell:
        "touch {output}"
    
