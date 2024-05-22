# swgs-absolutecn pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10040893.svg)](https://doi.org/10.5281/zenodo.10040893)

## Description

### *Summary*

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach.

### *Detailed description*

This pipeline implements a method for generating absolute copy number profiles from shallow whole genome sequencing data utilising a multi-stage fitting process to generate a set of read depth-normalised absolute copy number profiles. Copy number fits are generated in read count space using a modified implementation of [QDNAseq](https://github.com/ccagc/QDNAseq) on full read depth shallow whole genome BAM files. After which, a grid search of ploidy and purity values is performed to generate a matrix of potential absolute copy number profiles. The fitting of absolute copy number profiles includes an error function metric and a variant allele fraction anchoring process which assist in the selection of an optimal copy number fit. 

The error function, `clonality`, which measures segment residual from integer state across the genome and is computed as mean absolute error (MAE), in which a lower value is suggestive of better confirmation to integer copy number states, and therefore a better absolute copy number profile. variant allele fraction achoring utilises the near-ubiquitous presence of somatic homozygous _TP53_ variants in high grade serous ovarian cancer as an anchoring point for a copy number fit. Experimentally validated _TP53_ allele fractions at homozygous sites are compared to an expected homozygous _TP53_ variant allele fraction, as determined by the copy number value of the segment overlapping the _TP53_ locus. The smaller absolute differences between the experimental and expected _TP53_ variant allele fractions is suggestive of a better fit. Minimisation of the `clonality` error function, expected-to-experimental TP53 distance, and manual review of fits allows for the selection of an optimal absolute copy number profile for a given sample.

The selected fits are subject to a calculation of power to detect copy number alterations, as described [here](https://gmacintyre.shinyapps.io/sWGS_power/), in which the read depth of a given sample is assessed against its selected ploidy-purity combination. This process both determines if 1) a sample has a sufficient number of reads to support the selected ploidy-purity combination and 2) an optimal target number of reads to perform sample-specific read down sampling. This process normalises the read depth between samples in a ploidy and purity dependent manner so that read variance across segments is consistent, while excluding samples which are not supported by a sufficient number of reads.

Samples passing all filtering criteria then undergo read downsampling to the specified target number of reads determined by the previous steps and absolute copy number profiles fitted at the ploidy-purity combination selected prior.

## Table of contents

* [Pipeline setup](#pipeline-setup)
  + [Step 1 Clone the repo](#step-1-clone-the-repo)
  + [Step 2 Installing environment](#step-2-installing-environment)
  + [Step 3 Preparing the input files](#step-3-preparing-the-input-files)
    - [sample sheet](#sample-sheet)
    - [config.yaml](#configyaml)
    - [profile config.yaml](#profile-configyaml)
    - [workflow management](#workflow-management)
    - [updating the pipeline configuration](#updating-the-pipeline-configuration)
* [Running the pipeline](#running-the-pipeline)
  + [Step 4 Stage 1](#step-4-stage-1)
  + [Step 5 QC1](#step-5-qc1)
    - [Smoothing](#smoothing)
    - [Fit selection](#fit-selection)
  + [Step 6 Stage 2](#step-6-stage-2)
  + [Downsampling only](#downsampling-only)
* [Output files](#output-files)

## Compatibility

### Reference genome

This pipeline is currently only compatible with BAM files aligned with `hg19` or `GRCh37` genomes. There are no checks in place for using files aligned to an unsupported reference genome.

## Pipeline setup

### Step 1 Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
git clone https://github.com/Phil9S/swgs-absolutecn.git
cd swgs-absolutecn/
```

### Step 2 Installing environment

This pipeline can utilise either a conda environment or a containerised docker/singularity implementation to manage the software packages and dependecies.

- For the conda implementation please make sure either mamba or conda are installed and available on your system. Our recommendation is to use micromamba or, ideally, the installed conda version should utilise the libmamba solver library as the required environment contains a large number of packages. See installing [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) or [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
- For the docker/singularity implementation please make sure a compatible version of both snakemake and singularity are available on your system. See installing [singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

#### Conda environment
With conda or micromamba installed, from within the repository directory, run the `install_env.sh` script to generate a conda environment and install custom packages:
```
./install_env.sh mamba
```
or
```
./install_env.sh conda $HOME/miniconda/
```

If you are useing an existing conda installation, please use the conda installation directory relevant to your existing installation when running this section instead of `$HOME/miniconda/` to correctly initialise the conda environment.

The newly installed environment can be activated using the following:

```
micromamba activate swgs-abscn
```
or
```
conda activate swgs-abscn
```

#### Singularity-based implementation
If singularity and snakemake are already available, the pipeline can be run using a container by running the following:

Set the input and output directories (these should match those specified in the config.yaml - see step 3)
```
INPUTDIR="/path/to/inputfiles/"
OUTPUTDIR="/path/to/output/"
```
and then set the singularity args variable
```
SIGBIND="--use-singularity --singularity-args \"--bind ${INPUTDIR},${OUTPUTDIR}\""
```
The containerised environment can be found on docker hub here [phil9s/swgs-absolutecn](https://hub.docker.com/repository/docker/phil9s/swgs-absolutecn/)

### Step 3 Preparing the input files
#### Sample sheet

The workflow requires a single input file `sample_sheet.tsv` which is a tab-separated document detailing sample names, patient names, smoothing booleans, TP53 allele frequencies, and BAM file locations. Make sure that the BAM file paths are absolute and that sample and patient names are sensible (i.e. do not contain abstract characters or white space). The `TPfreq` column should only contain a float (range 0.00-1.00) or `NA`. The `smooth` column is boolean and should only contain either `TRUE` or `FALSE`. Optionally, users may provide precomputed or orthogonally generated sample ploidy and purity estimates using the `precPloidy` and `precPurity` fields but these are not required. Both `precPloidy` and `prePurity` can be provided seperately or together and samples do not need to have either field consistently, meaning one sample may have ploidy/purity estimate and others may not and the pipeline will only perform grid searches across the missing values (e.g In the example table below, a full gridsearch will occur for SAM1, no search performed for SAM2, only ploidy for SAM3, and only purity for SAM4). 

The table below demonstrates the basic schema required and the workflow will validate this file prior to running.

|PATIENT_ID|SAMPLE_ID|TP53freq|smooth|file          |precPloidy|precPurity|
|----------|---------|--------|------|--------------|----------|----------|
|PAT1      |SAM1     |0.45    |TRUE  |/data/SAM1.bam|NA        |NA        |
|PAT1      |SAM2     |0.55    |FALSE |/data/SAM2.bam|2.4       |0.75      |
|PAT2      |SAM3     |NA      |FALSE |/data/SAM3.bam|NA        |0.43      |
|PAT2      |SAM4     |NA      |TRUE  |/data/SAM4.bam|3.2       |NA        |

An example sample_sheet.tsv is included in this repository.

#### config.yaml

The config.yaml (`config/config.yaml`) contains the necessary information for the pipeline you wish to run. This includes the location of the samplesheet.tsv, bin size, project name, output directory, and filtering parameters.

#### Filtering parameters

Various filters for acceptable ploidy/purity combinations for fitting absolute copy number can be modified or disabled depending user requirements.

|variable            |default |function                                                                                                                 |type       |values        |
|--------------------|--------|-------------------------------------------------------------------------------------------------------------------------|-----------|--------------|
|af_cutoff           |0.15    |Maximum difference between expected and observed _TP53_ allele fraction                                                  |float      |0.0-1.0       |
|use_seed            |"TRUE"  |Set whether to use `seed_val` to ensure CBS segmentation returns identical results                                       |string bool|"TRUE","FALSE"|
|seed_val            |"9999"  |seed value used by `use_seed`                                                                                            |string     |any string    |
|filter_underpowered |"TRUE"  |Set whether to filter ploidy/purity combinations without sufficent available read depth to support the given profile     |string bool|"TRUE","FALSE"|
|ploidy_min          |1.6     |Minimum ploidy value for gridsearch - must be less than `ploidy_max`. Ignored for a sample if `precPloidy` is provided   |float      |1.0-20.0      |
|ploidy_max          |8.0     |Minimum ploidy value for gridsearch - must be greater than `ploidy_min`. Ignored for a sample if `precPloidy` is provided|float      |1.0-20.0      |
|purity_min          |0.15    |Minimum purity value for gridsearch - must be less than `purity_max`. Ignored for a sample if `precPurity` is provided   |float      |0.0-1.0       |
|purity_max          |1.0     |Minimum purity value for gridsearch - must be greater than `purity_min`. Ignored for a sample if `precPurity` is provided|float      |0.0-1.0       |
|filter_homozygous   |"TRUE"  |Set whether to filter ploidy/purity combinations with a proportion of homozygous loss greater than `homozygous_prop`     |string bool|"TRUE","FALSE |
|homozygous_prop     |10000000|Proportion of genome (in basepairs) at which to filter a ploidy/purity combination where `filter_homozygous` is "TRUE"   |integer    |minimum=0     |
|homozygous_threshold|0.4     |Threshold at which to assign homozygous loss to copy number segment counted by `homozygous_prop`                         |float      |0.0-0.99      |

For most users, the default parameters should work well but in certain instances, these values should be modifed. 

For example, if users are not generating a sufficent number of high quality absolute copy number profiles, setting `filter_underpowered` to "FALSE" will show a larger range of fits which, although statistically underpowered, could be correct for given sample.

Another use case would be samples where the general purity range is known, for example cell line or LCM data, where the expected purity is known to be ~1.0. As such setting `purity_min` to 0.95 would restrict ploidy/purity combinations to enforce this expectation.

#### Workflow management
##### profile config.yaml

The profile configs (cluster_config.yaml & config.yaml) (`profile/*/*`) contains the necessary information to configure the job submission parameters passed to a given workload manager (or lack thereof). This includes the number of concurrently sumbitted jobs, account/project name, partition/queue name, and default job resources (though these are low and should work on almost any cluster). This pipeline was primarily developed using the [SLURM](https://slurm.schedmd.com/documentation.html) work load manager for job submission by snakemake. For individuals running on non-workload managed clusters, or utilising other workload managers, profiles are provided to allow for job submission with minimal configuration. Currently supported profiles are `local`, `slurm`, and `pbs`. These can be edited via the script described in the next section.

#### Updating the pipeline configuration

Environment-specific and pipeline-specific parameters need to be set for each run of this pipeline. While it is possible to manually edit the YAML files, a script has been provided to update the most frequently altered parameters programmatically. The script will iteratively list the parameter and its current value, asking for a user submitted new value should it be needed. If the value is already acceptable or does not need to be changed then an empty value (enter return without typing) will keep the current setting.

Run one of the following code for the profile you wish to update (typically both the pipeline configuration and one cluster configuration):

```
# conda - with swgs-abscn env
# Pipeline configuration
./update_configs.py -c config
# Pipeline filters
./update_configs.py -c filters
# workflow configuration
./update_configs.py -c {slurm,pbs,local}
```
singularity users can run the following (whilst within the repository directory):
```
SIGCMD="singularity exec --bind "$(pwd -P)" docker://phil9s/swgs-absolutecn:latest $(pwd -P)"
# Pipeline configuration
$SIGCMD/update_configs.py -c config
# Pipeline filters
$SIGCMD/update_configs.py -c filters
# workflow configuration
$SIGCMD/update_configs.py -c {slurm,pbs,local}

```

## Running the pipeline

### Step 4 Stage 1

Once the pipeline and cluster parameters have been set and the samplesheet is prepared, the first stage of the pipeline is ready to run.
Be sure to update the profile path to match your cluster/server configuration. 

Confirm the pipeline is configured correctly, run the first stage using the `dry-run` mode.

```
# conda - with swgs-abscn env
snakemake -n --profile profile/slurm/ --snakefile stage_1
```
```
# singularity
snakemake -n --profile profile/slurm/ --snakefile stage_1 ${SIGBIND}
```
If the previous step ran without error then run the following:
```
# conda - with swgs-abscn env
snakemake --profile profile/slurm/ --snakefile stage_1
```
```
# singularity
snakemake --profile profile/slurm/ --snakefile stage_1 ${SIGBIND}
```

### Step 5 QC1

At the conclusion of stage 1, files and fits will be generated for all samples present in the `samplesheet.tsv` provided. This will include grid search-generated fits, copy number profile plots, and a `QDNAseq` RDS file containing the copy number fit data. This data is not immediately usable in downstream processes and fits must undergo quality control and fit selection, as samples may generate more than one viable copy number profile. 

#### Smoothing

Prior to fit selection, a subset of samples may require smoothing of segments in order to be viable. Read the guide provided here to select and update which samples require smoothing [here](resources/smoothing_guide.md). Once the `samplesheet.tsv` has been updated with the new `smooth` values, rerun stage 1. Be sure to update the profile path to match your cluster/server configuration.

```
# conda - with swgs-abscn env
snakemake --profile profile/slurm/ -F all --snakefile stage_1
```
```
#singularity
snakemake --profile profile/slurm/ -F all --snakefile stage_1 ${SIGBIND}
```

#### Fit selection

After running stage 1 with the appropraite smoothing values, copy number fits will have been generated for each sample. In most cases, multiple viable fits will have been selected for each sample and a semi-qualitative quality control process needs to be applied to select the best fitting solution or exclude a sample should no fit be good enough.

Follow the guide on fit selection [here](resources/quality_control_guide.md) to perform quality control and update the `{project}_fit_QC_predownsample.tsv` file. Once this step has been performed and a single acceptable fit has been selected for each sample, proceed to stage 2.

### Step 6 Stage 2

Provided quality control and fit selection was performed correctly, stage 2 of the pipeline can be performed which refits all copy number profiles using downsampled BAM files and the selected fits from stage 1.
Be sure to update the profile path to match your cluster/server configuration. 

As before, confirm the pipeline is configured correctly by running with the `dry-run` mode.

```
# conda - with swgs-abscn env
snakemake -n --profile profile/slurm/ --snakefile stage_2
```
```
# singularity
snakemake -n --profile profile/slurm/ --snakefile stage_2 ${SIGBIND}
```
and if the previous step ran without error then run the following:
```
# conda - with swgs-abscn env
snakemake --profile profile/slurm/ --snakefile stage_2
```
```
# singularity
snakemake --profile profile/slurm/ --snakefile stage_2 ${SIGBIND}
```

### Downsampling only

In some instances, where ploidy and purity are provided from alternative sources, users may only want to perform the read depth normalisation step and not perform a gridsearch. In these cases, if all samples in the `sample_sheet.tsv` contain a valid `precPloidy` and `precPurity` values the `processPrecomputed.R` script can generate all the required output files to start stage 2 without needing to perform any gridsearch or select a ploidy/purity combination, and instead skip directly to read depth normalisation (via read downsampling) using the provied `precPloidy` and `precPurity` combination.

First, check all samples have valid `precPloidy` and `precPurity` values by running the following:
```
# conda - with swgs-abscn env
snakemake -n --profile profile/slurm/ --snakefile stage_1
```
```
# singularity
snakemake -n --profile profile/slurm/ --snakefile stage_1 ${SIGBIND}
```
This should return an error message stating `all ploidy and purity values are precomputed - Run Rscript scripts/processPrecomputed.R and skip to stage_2`. If this did not occur and stage 1 ran a dry-run then not all samples have a `precPloidy` and `precPurity` value specified.

Provided that the previous step validated and returned the specified message then the following can be run to generate the required output files to skip stage 1

```
# conda - with swgs-abscn env
Rscript scripts/processPrecomputed.R
```
```
# singularity
SIGCMD="singularity exec --bind "$(pwd -P)" --bind ${INPUTDIR},${OUTPUTDIR} docker://phil9s/swgs-absolutecn:latest"
$SIGCMD Rscript scripts/processPrecomputed.R
```

At which point Stage 2 can then be performed as previously described.

## Output files

The final output is generated by stage 2 is three files located in the specified output directory `{output_dir}/sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/`
- `*_ds_abs_fits.tsv` - Tab-seperated file containing absolute copy number profile metadata and fitting information
- `*_ds_absCopyNumber.rds` - QDNASeqmod class object containing binned copy number data in `rds` format
- `*_ds_absCopyNumber_segTable.tsv` - A multisample segment table consisting of copy number segment coordinates and associated segment values as unrounded absolute copy number

## Authors

* Philip Smith (@phil9s)

## Citation

Please cite this pipeline using the publication ["The copy-number landscape of recurrent ovarian high grade serous carcinoma"](https://doi.org/10.1038/s41467-023-39867-7) Smith & Bradley et al. 2023 - _Nature Communications_ and/or the version controlled zenodo repository [10.5281/zenodo.10040893](https://zenodo.org/doi/10.5281/zenodo.10040893).

