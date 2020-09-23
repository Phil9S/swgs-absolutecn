## Smoothing guide

In most cases, no prior information on whether the boolean `smooth` column in the `samplesheet.tsv` should contain a `TRUE` or `FALSE` for each sample and so by default should be set to `FALSE`. For most samples, smoothing is not required but oversegmentation issues and data-specific factors mean that a subset of samples will likely require smoothing to be performed and as such `smooth` to be set to `TRUE`.

To perform an assessment of which samples require smoothing, plots generated in stage 1 can be visually assessed. The files are outuput into the directory `{output_dir}/sWGS_fitting/{project}_{bin}kb/absolute_PRE_down_sampling/plots/` where `output_dir`, `project`, and `bin` corespond to the values set when the pipeline was run. (DOUBLE CHECK DIR PATH)

Additionally, assessing the total number of called segments (annotated above each copy number profile plot) can be a reasonable indicator as well, though this will be tissue/cancer specific. To provide some context for which samples do or do not need smoothing, the images below are examples of ovarian cancer copy number profiles which would benefit from having smoothing applied.

![Smoothing example 1](resources/images/smoothing_example_1.png)
![Smoothing example 2](resources/images/smoothing_example_2.png)

Samples with copy number fits resembling the fits shown above should most likely require smoothing. After updating the `smooth` boolean for samples requiring smoothing, rerun stage 1 using the newly edited `samplesheet.tsv`, as described below:
