## added conditionals for CRAM support
if config["filetype"] in ["CRAM"]:

    rule cram_to_bam:
        input:
            bam=get_bam,
            check=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"),
        output:
            bam=temp(Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam")),
            bai=temp(Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam.bai")),
        log:
            "logs/cram_to_bam_{project}_{bin}kb_{sample}.log",
        container:
            image
        params:
            reference=config["reference"],
        shell:
            "(samtools view -b -T {params.reference} {input.bam} > {output.bam} && samtools index {output.bam} {output.bai}) 2> {log}"

else:
    # symlinking run locally
    localrules:
        sym_link,

    rule sym_link:
        input:
            bam=get_bam,
            check=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"),
        output:
            Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam"),
        params:
            bam_abs=lambda wc, input: os.path.abspath(input.bam),
        log:
            "logs/sym_link_{project}_{bin}kb_{sample}.log",
        container:
            image
        shell:
            "ln -s {params.bam_abs} {output} 2> {log}"
