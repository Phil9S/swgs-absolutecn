## added conditionals for CRAM support
if config["filetype"] in ["CRAM"]:

    rule cram_to_bam:
        input:
            bam=get_bam,
            check=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"),
        output:
            temp(
                expand(
                    Path(OUT_DIR,"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
                )
            ),
        log:
            "logs/cram_to_bam_{project}_{bin}kb_{sample}.log",
        container:
            image
        threads: 1
        params:
            reference=config["reference"],
        shell:
            "samtools view -b -T {params.reference} {input.bam} > {output} && samtools index {output}"

else:

    rule sym_link:
        input:
            bam=get_bam,
            check=Path(OUT_DIR,"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"),
        output:
            expand(Path(OUT_DIR,"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")),
        log:
            "logs/sym_link_{project}_{bin}kb_{sample}.log",
        container:
            image
        threads: 1
        shell:
            "ln -s {input.bam} {output}"
