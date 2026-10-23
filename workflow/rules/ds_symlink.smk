# symlinking run locally
localrules:
    ds_symlink,

rule ds_symlink:
    input:
        bam=get_bam,
    output:
        Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam"),
    params:
        bam_abs=lambda wc, input: os.path.abspath(input.bam),
    log:
        "logs/ds_symlink_{project}_{bin}kb_{sample}.log",
    container:
        image
    shell:
        "ln -s {params.bam_abs} {output} 2> {log}"
