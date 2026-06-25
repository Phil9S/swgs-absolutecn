## added conditionals for CRAM support
if config['filetype'] in ["CRAM"]:
    rule cram_to_bam:
        input:
            bam=get_bam, 
            check=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"
        params:
            reference=config["reference"]
        output:
            temp(expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam"))
        log:
            "logs/cram_to_bam_{project}_{bin}kb_{sample}.log"
        threads: 1
        container:
            image
        shell:
            "samtools view -b -T {params.reference} {input.bam} > {output} && samtools index {output}"
else:
    rule sym_link:
        input:
            bam=get_bam,
            check=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}_bam.ok"
        output:
            expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
        log:
            "logs/sym_link_{project}_{bin}kb_{sample}.log"
        threads: 1
        container:
            image
        shell:
            "ln -s {input.bam} {output}"
