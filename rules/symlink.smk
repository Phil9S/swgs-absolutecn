rule sym_link:
    input:
        bam=get_bam,
        check=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bam.ok"
    output:
        expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
    threads: 1
    shell:
        "ln -s {input.bam} {output}"
