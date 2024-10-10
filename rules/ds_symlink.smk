rule ds_symlink:
        input:
            bam=get_bam,
        output:
            expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
        threads: 1
        shell:
            "ln -s {input.bam} {output}"

