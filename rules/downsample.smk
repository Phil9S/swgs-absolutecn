rule downsample:
    input:
        bam=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/bams/{sample}.bam",
        meta=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_fit_QC_predownsample.tsv",
        rds=OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/{project}_{bin}kb_relSmoothedCN.rds"
    output:
        OUT_DIR+"sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/downsampled_bams/{sample}.bam"
    threads: 1
    params:
        outdir=OUT_DIR,
        bin="{bin}",
        project="{project}"
    shell:
        "Rscript scripts/downsampleBams.R {input.bam} {wildcards.sample} {input.meta} {input.rds} {params.outdir} {params.bin} {params.project}"
    
