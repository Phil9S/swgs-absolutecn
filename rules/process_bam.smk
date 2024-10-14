## added conditionals for CRAM support
if config['filetype'] in ["CRAM"]:
    rule cram_to_bam:
        input:
            bam=get_bam,
            check=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bam.ok"
        params:
            reference=config["reference"]
        output:
            temp(expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam"))
        threads: 1
        singularity:
            image_base_url+"swgs-absolutecn:latest"
        shell:
            "samtools view -b -T {params.reference} {input.bam} > {output} && samtools index {output}"
else:
    rule sym_link:
        input:
            bam=get_bam,
            check=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bam.ok"
        output:
            expand(OUT_DIR+"sWGS_fitting/{{project}}_{{bin}}kb/bams/{{sample}}.bam")
        threads: 1
        shell:
            "ln -s {input.bam} {output}"
