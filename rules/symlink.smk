rule sym_link:
    input:
        bam=get_bam
    output:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/full_rd/bams/{{sample}}.bam")
        #expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/full_rd/bams/{{sample}}.bai")
    params:
        bin=config["bins"]
    shell:
        "#ls {input}"
        "#touch {output}"
