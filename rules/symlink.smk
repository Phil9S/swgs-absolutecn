rule sym_link:
    input:
        bam=get_bam
    output:
        expand(OUT_DIR+"sWGS_fitting/bams/{{sample}}.bam")
    params:
        bin=config["bins"]
    shell:
        "ln -s {input} {output}"
