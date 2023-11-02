# swgs-absolutecn pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10040893.svg)](https://doi.org/10.5281/zenodo.10040893)

## Description

*Summary*

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach.

*Detailed description*

This pipeline implements a method for generating absolute copy number profiles from shallow whole genome sequencing data utilising a multi-stage fitting process to generate a set of read depth-normalised absolute copy number profiles. Copy number fits are generated in read count space using a modified implementation of [QDNAseq](https://github.com/ccagc/QDNAseq) on full read depth shallow whole genome BAM files. After which, a grid search of ploidy and purity values is performed to generate a matrix of potential absolute copy number profiles. The fitting of absolute copy number profiles includes an error function metric and a variant allele fraction anchoring process which assist in the selection of an optimal copy number fit. 

The error function, `clonality`, which measures segment residual from integer state across the genome and acts similarly to a function such as RMSE, in which a lower value is suggestive of better confirmation to integer copy number states, and therefore a better absolute copy number profile. variant allele fraction achoring utilises the near-ubiquitous presence of somatic homozygous _TP53_ variants in high grade serous ovarian cancer as an anchoring point for a copy number fit. Experimentally validated _TP53_ allele fractions at homozygous sites are compared to an expected homozygous _TP53_ variant allele fraction, as determined by the copy number value of the segment overlapping the _TP53_ locus. The smaller absolute differences between the experimental and expected _TP53_ variant allele fractions is suggestive of a better fit. Minimisation of the `clonality` error function, expected-to-experimental TP53 distance, and manual review of fits allows for the selection of an optimal absolute copy number profile for a given sample.

The selected fits are subject to a calculation of power to detect copy number alterations, as described [here](https://gmacintyre.shinyapps.io/sWGS_power/), in which the read depth of a given sample is assessed against its selected ploidy-purity combination. This process both determines if 1) a sample has a sufficient number of reads to support the selected ploidy-purity combination and 2) an optimal target number of reads to perform sample-specific read down sampling. This process normalises the read depth between samples in a ploidy and purity dependent manner so that read variance across segments is consistent, while excluding samples which are not supported by a sufficient number of reads.

Samples passing all filtering criteria then undergo read downsampling to the specified target number of reads determined by the previous steps and absolute copy number profiles fitted at the ploidy-purity combination selected prior.

## Table of contents

* [Pipeline setup](#pipeline-setup)
  + [Step 1 Clone the repo](#step-1-clone-the-repo)
  + [Step 2 Install conda](#step-2-install-conda)
  + [Step 3 Installing environment](#step-3-installing-environment)
  + [Step 4 Preparing the input files](#step-4-preparing-the-input-files)
    - [Sample sheet](#sample-sheet)
    - [config.yaml](#configyaml)
    - [profile config.yaml](#profile-configyaml)
    - [workflow management](#workflow-management)
    - [Updating the pipeline configuration](#updating-the-pipeline-configuration)
* [Running the pipeline](#running-the-pipeline)
  + [workflow management](#workflow-management)
  + [Step 5 Stage 1](#step-5-stage-1)
  + [Step 6 QC1](#step-6-qc1)
    - [Smoothing](#smoothing)
    - [Fit selection](#fit-selection)
  + [Step 7 Stage 2](#step-7-stage-2)
  + [Step 8 QC2](#step-8-qc2)

## Compatibility

### Reference genome

This pipeline is currently only compatible with BAM files aligned with `hg19` or `GRCh37` genomes.
There are no checks in place for using the correct genome and results will be likely be extremely poor or incorrect if done using an unsupported reference genome.

## Pipeline setup

### Step 1 Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
git clone https://github.com/Phil9S/swgs-absolutecn.git
cd swgs-absolutecn/
```

### Step 2 Install conda

This pipeline utilises micromamba or conda installation environments to manage the software packages and versions. Please make sure either mamba or conda are installed and available on your system.
Our recommendation is to use micromamba or, ideally, the installed conda version should utilise the libmamba solver library as the required environment contains a large number of packages.

See installing [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) or [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for more information.

### Step 3 Installing environment

From within the repository directory, run the `install_env.sh` script to generate a conda environment and install custom packages:

```
./install_env.sh mamba
```
or
```
install_env.sh conda $HOME/miniconda/
```

If you used a previously installed conda build please use the conda or miniconda installation directory when running this section instead of '$HOME/miniconda/' to correctly initialise the conda environment.

The newly installed environment can be activated using the following:

```
micromamba activate swgs-abscn
```
or
```
conda activate swgs-abscn
```

### Step 4 Preparing the input files

#### Sample sheet

The workflow requires a single input file `sample_sheet.tsv` which is a tab-separated document detailing sample names, patient names, smoothing booleans, TP53 allele frequencies, and BAM file locations. Make sure that the BAM file paths are absolute and that sample and patient names are sensible (i.e. do not contain abstract characters or white space). The `TPfreq` column should only contain a float (range 0.00-1.00) or `NA`. The `smooth` column is boolean and should only contain either `TRUE` or `FALSE`.

The table below demonstrates the basic schema required and the workflow will validate this file prior to running.

|PATIENT_ID|SAMPLE_ID|TP53freq|smooth|file         |
|----------|---------|--------|------|-------------|
|PAT1      |SAM1     |0.45    |TRUE  |/data/SAM1.bam|
|PAT1      |SAM2     |0.55    |FALSE |/data/SAM2.bam|

An example sample_sheet.tsv is included in this repository.

#### config.yaml

The config.yaml (`config/config.yaml`) contains the necessary information for the pipeline you wish to run. This includes the location of the samplesheet.tsv, bin size, project name, and output directory.

#### profile config.yaml

The profile configs (cluster_config.yaml & config.yaml) (`profile/*/*`) contains the necessary information to configure the job submission parameters passed to a given workload manager (or lack thereof). This includes the number of concurrently sumbitted jobs, account/project name, partition/queue name, and default job resources (though these are low and should work on almost any cluster).

#### Workflow management

This pipeline was primarily developed using the [SLURM](https://slurm.schedmd.com/documentation.html) work load manager for job submission by snakemake. For individuals running on non-workload managed clusters, or utilising other workload managers, profiles are provided to allow for job submission with minimal configuration. Currently supported profiles are `local`, `slurm`, and `pbs`. These can be edited via the script described in the next section.

#### Updating the pipeline configuration

Environment-specific and pipeline-specific parameters need to be set for each run of this pipeline. While it is possible to manually edit the YAML files, a script has been provided to update the most frequently altered parameters programmatically. The script will iteratively list the parameter and its current value, asking for a user submitted new value should it be needed. If the value is already acceptable or does not need to be changed then an empty value (enter return without typing) will keep the current setting.

With the `swgs-abscn` conda environment active, run one of the following code for the profile you wish to update (typically both the pipeline configuration and one cluster configuration):

```
# Pipeline configuration
python update_configs.py -c config
# Slurm configuration
python update_configs.py -c slurm
# PBS-torque configuration
python update_configs.py -c pbs
# Local/non-managed configuration
python update_configs.py -c local
```

## Running the pipeline

### Step 5 Stage 1

Once the pipeline and cluster parameters have been set and the samplesheet is prepared, the first stage of the pipeline is ready to run.

With the `swgs-abscn` conda environment active, run the following:

Confirm the pipeline is configured correctly, run the first stage using the `dry-run` mode.

```
snakemake -n --profile profile/*/ --snakefile stage_1
```

If the previous step ran without error then run the following:

```
snakemake --profile profile/*/ --snakefile stage_1
```

Where * is replaced with the profile matching your cluster/server configuration. 

### Step 6 QC1

At the conclusion of stage 1, files and fits will be generated for all samples present in the `samplesheet.tsv` provided. This will include grid search-generated fits, copy number profile plots, and a `QDNAseq` RDS file containing the copy number fit data. This data is not immediately usable in downstream processes and fits must undergo quality control and fit selection, as samples may generate more than one viable copy number profile.

#### Smoothing

Prior to fit selection, a subset of samples will likely require smoothing of segments in order to be viable. Read the guide provided here to select and update which samples require smoothing [here](resources/smoothing_guide.md). Once the `samplesheet.tsv` has been updated with the new `smooth` values, rerun stage 1 using the following:

```
snakemake --profile profile/*/ -F all --snakefile stage_1
```

Where * is replaced with the profile matching your cluster/server configuration.

#### Fit selection

After running stage 1 with the appropraite smoothing values, copy number fits will have been generated for each sample. In most cases, multiple viable fits will have been selected for each sample and a semi-qualitative quality control process needs to be applied to select the best fitting solution or exclude a sample should no fit be good enough.

Follow the guide on fit selection [here](resources/quality_control_guide.md) to perform quality control and update the `{project}__fit_QC_predownsample.tsv` file. Once this step has been performed and a single acceptable fit has been selected for each sample, proceed to step 7.

### Step 7 Stage 2

Provided quality control and fit selection was performed correctly, stage 2 of the pipeline can be performed which refits all copy number profiles using downsampled BAM files and the selected fits from stage 1.

As before, confirm the pipeline is configured correctly by running with the `dry-run` mode.

```
snakemake -n --profile profile/*/ --snakefile stage_2
```

and if the previous step ran without error then run the following:

```
snakemake --profile profile/*/ --snakefile stage_2
```

Where * is replaced with the profile matching your cluster/server configuration.

### Step 8 QC2

To confirm the quality of newly generated downsampled absolute copy number profiles generated by stage 2, an evalution of outputted fits should be performed as described previously in step 6 [here](resources/quality_control_guide.md). The output `{}` should be updated accordingly, with poor fits being excluded.

## Authors

* Philip Smith (@phil9s)

## Citation

Please cite this pipeline using the publication ["The copy-number landscape of recurrent ovarian high grade serous carcinoma"](https://doi.org/10.1038/s41467-023-39867-7) Smith & Bradley et al. 2023 - _Nature Communications_ and/or the version controlled zenodo repository [10.5281/zenodo.10040893](https://zenodo.org/doi/10.5281/zenodo.10040893).

