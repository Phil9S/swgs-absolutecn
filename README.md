# swgs-absolutecn pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10040893.svg)](https://doi.org/10.5281/zenodo.10040893)

## Description

### *Summary*

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach using a snakemake-based workflow.

<details>
<summary>Detailed description</summary>

### *Detailed description*

This pipeline implements a method for generating absolute copy number profiles from shallow whole genome sequencing data utilising a multi-stage fitting process to generate a set of read depth-normalised absolute copy number profiles. Copy number fits are generated in read count space using a modified implementation of [QDNAseq](https://github.com/ccagc/QDNAseq) on full read depth shallow whole genome BAM files. After which, a grid search of ploidy and purity values is performed to generate a matrix of potential absolute copy number profiles. The fitting of absolute copy number profiles includes an error function metric and a variant allele fraction anchoring process which assist in the selection of an optimal copy number fit. 

The error function, `clonality`, which measures segment residual from integer state across the genome and is computed as mean absolute error (MAE), in which a lower value is suggestive of better confirmation to integer copy number states, and therefore a better absolute copy number profile. variant allele fraction achoring utilises the near-ubiquitous presence of somatic homozygous _TP53_ variants in high grade serous ovarian cancer as an anchoring point for a copy number fit. Experimentally validated _TP53_ allele fractions at homozygous sites are compared to an expected homozygous _TP53_ variant allele fraction, as determined by the copy number value of the segment overlapping the _TP53_ locus. The smaller absolute differences between the experimental and expected _TP53_ variant allele fractions is suggestive of a better fit. Minimisation of the `clonality` error function, expected-to-experimental TP53 distance, and manual review of fits allows for the selection of an optimal absolute copy number profile for a given sample.

The selected fits are subject to a calculation of power to detect copy number alterations, as described [here](https://gmacintyre.shinyapps.io/sWGS_power/), in which the read depth of a given sample is assessed against its selected ploidy-purity combination. This process both determines if 1) a sample has a sufficient number of reads to support the selected ploidy-purity combination and 2) an optimal target number of reads to perform sample-specific read down sampling. This process normalises the read depth between samples in a ploidy and purity dependent manner so that read variance across segments is consistent, while excluding samples which are not supported by a sufficient number of reads.

Samples passing all filtering criteria then undergo read downsampling to the specified target number of reads determined by the previous steps and absolute copy number profiles fitted at the ploidy-purity combination selected prior.
</details>

## Pipeline setup

### Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your  system.

```
git clone https://github.com/Phil9S/swgs-absolutecn.git
cd swgs-absolutecn/
```
### Install environment

This pipeline can utilise either a virtual environment or a container to manage the software packages and dependecies.

#### Package manager

<details>
<summary>pixi (recommended)</summary>

Follow the linked instructions to install [Pixi](https://pixi.prefix.dev/latest/installation/) then run the following:

```
pixi install --run-post-link-scripts && pixi run ./install_env.sh
```
</details>

<details>
<summary>micromamba</summary>

Follow the linked instructions to install [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/) then run the following:

```
micromamba create -f conda.yaml
micromamba activate swgs-abscn && ./install_env.sh
```
</details>

<details>
<summary>conda</summary>

Follow the linked instructions to install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) then run the following:

```
conda create -f conda.yaml
conda activate swgs-abscn && ./install_env.sh
```
</details>

#### Container-based

For the containerised implementation, a compatible version of both snakemake and singularity are required on your system. See installing [singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). This can done using one of the previous methods (pixi, micromamba, conda, etc.), included with the system via modules or, installed directly.

Prior to using the containerised implementation, set the following `ENV` variables to bind the appropriate volumes.

Set the input and output directories (these should match those specified in the config.yaml)
```
INPUTDIR="/path/to/inputfiles/"
OUTPUTDIR="/path/to/output/"
```
set the singularity args variable
```
SIGBIND="--use-singularity --singularity-args \"--bind ${INPUTDIR},${OUTPUTDIR}\""
```
as well as the cmdline
```
SIGCMD="singularity exec --bind "$(pwd -P)" --bind ${INPUTDIR},${OUTPUTDIR} docker://phil9s/swgs-absolutecn:latest"
```
The container can be found on docker hub [here](https://hub.docker.com/r/phil9s/swgs-absolutecn).

### Preparing workflow
#### Sample sheet

The workflow requires a single tab-separated file, `sample_sheet.tsv`, which specifies the samples and associated BAM files to use as inputs for the workflow. The `sample_sheet.tsv` schema is validated prior to executing the workflow.

- BAM file location, `file`, should be an absolute path.
- `SAMPLE_ID` must be unique.
- `TPfreq` column is an optional allele frequency for _TP53_ (float (0-1.00) or `NA`).
- The `smooth` column is boolean (`TRUE` or `FALSE`, also see `config.yaml` for `autosmoothing`). 
- Users may provide precomputed ploidy and purity estimates, `precPloidy` and `precPurity` fields and both values are not required per sample.

##### Sample sheet schema

|PATIENT_ID|SAMPLE_ID|TP53freq|smooth|file          |precPloidy|precPurity|
|----------|---------|--------|------|--------------|----------|----------|
|PAT1      |SAM1     |0.45    |TRUE  |/data/SAM1.bam|NA        |NA        |
|PAT1      |SAM2     |0.55    |FALSE |/data/SAM2.bam|2.4       |0.75      |
|PAT2      |SAM3     |NA      |FALSE |/data/SAM3.bam|NA        |0.43      |
|PAT2      |SAM4     |NA      |TRUE  |/data/SAM4.bam|3.2       |NA        |

#### config.yaml

The config.yaml (`config/config.yaml`) contains the necessary information for the pipeline you wish to run and filters to apply. This file should be edited to match your desired configuration.

##### **fixed parameters**

|variable|default|function|type|values|
|--------|-------|--------|----|------|
|samplesheet|"sample_sheet.tsv"|`sample_sheet.tsv` location|string|file path|
|out_dir|"results/"|Root location for workflow output|string|directory path|
|project_name|"example_run"|Name for workflow run to include in outputs|string|any string|
|bins|30|Size of genomic bins (in kilobases)|integer|1,5,15,30,50,100,500,1000|
|use_seed|"TRUE"|Set whether to use `seed_val` to ensure CBS segmentation returns identical results|bool|"TRUE","FALSE"|
|seed_val|"9999"|Value of seed to use|string|any string|
|genome|"hg19"|Name of reference genome to use for QDNAseq bin annotations|string|"hg19","hg38" |
|filetype|"BAM"|File type of input file - either a aligned BAM or CRAM file|string|"BAM","CRAM"|
|reference|""|File path to the reference used in generation of CRAM files if using the filetype = "CRAM" option|string|file path|
|pairedEnd|"TRUE"|Boolean whether input files are paired or single end (only affects expected variance)|bool|"TRUE","FALSE"|

##### **automation parameters**

|variable|default|function|type|values|
|--------|-------|--------|----|------|
|autoSmooth|"FALSE"|Apply automatic smoothing of over-segmented copy number profiles|bool|"TRUE","FALSE"|
|smoothThreshold|700|Amount of segments to trigger auto smoothing|integer|44-9999|
|flagThreshold|0.84|auto fitting prediction probability threshold at which to flag samples as potentially poor|float|0.1.0|
fitMethod|"errorOnly"|Autofitting method to use to select best copy number profile fit from gridsearch|string|"errorOnly","randforest"|
|errorMetric|"clonality"|Error metric used in autofitting to determine `errorOnly` best fit or tiebreaking procedure in `randforest` autofitting methods|string|"clonality","segvariance","rmse"|

##### **filtering parameters**

Various filters for acceptable ploidy/purity combinations for fitting absolute copy number can be modified or disabled depending user requirements.

|variable|default|function|type|values|
|--------|-------|--------|----|------|
|af_cutoff|0.15|Maximum difference between expected and observed _TP53_ allele fraction|float|0.0-1.0|
|filter_underpowered |"TRUE"|Set whether to filter ploidy/purity combinations without sufficent available read depth to support the given profile|string bool|"TRUE","FALSE"|
|ploidy_min|1.6|Minimum ploidy value for gridsearch|float|1.0-20.0|
|ploidy_max|8.0|Minimum ploidy value for gridsearch|float|1.0-20.0|
|purity_min|0.15|Minimum purity value for gridsearch|float|0.0-1.0|
|purity_max|1.0|Minimum purity value for gridsearch|float|0.0-1.0|
|filter_homozygous|"TRUE"|Set whether to filter ploidy/purity combinations with a proportion of homozygous loss greater than `homozygous_prop`|bool|"TRUE","FALSE|
|homozygous_prop|10000000|Proportion of genome (in basepairs) at which to filter a ploidy/purity combination where `filter_homozygous` is "TRUE"|integer|minimum=0|
|homozygous_threshold|0.4|Threshold at which to assign homozygous loss to copy number segment counted by `homozygous_prop`|float|0.0-0.99|

#### Workflow management
##### profile config.yaml

The snakemake workflow is designed to run using job schedulers and as such uses the [snakemake executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/) to manage job submission and request compute resources. By default, the provided environments and containers include profiles for `slurm` and `lsf`, as well as  `local` implementation.

The profile configs (`profile/*/profile.yaml`) contains the necessary information to configure the job submission parameters passed to a given workload manager (or lack thereof). Edit the appropriate `profile.yaml` to match your given compute environment or add your own executor plugin as required.

## Running the pipeline

### Workflow validation

Validate the environment and pipeline is configured correctly by run the following using the `dry-run` mode (append the `${SIGBIND}` ENV variable if using a containerised implementation).

```
snakemake -n -p --profile profile/slurm/ --snakefile auto
```

### Runing workflow

The `swgs-absolutecn` workflow provides two approaches to copy number profile fitting, `Automated` and `Staged`.
 
- `Automated` runs the entire workflow, including segmentation, ploidy/purity gridsearch, absolute copy number fitting, read depth normalisation downsampling, utilising an autofitting procedure (error minimisation or random forest prediction) to select the best copy number fit.

- `Staged` splits the workflow into stage 1 and stage 2, where manual quality control and selection of copy number fit is performed by the user after the initial gridsearch is performed.

#### Automated workflow

<details>
<summary></summary>

```
# conda - with swgs-abscn env
snakemake --profile profile/slurm/ --snakefile auto
```

</details>

#### Staged workflow

<details>
<summary></summary>

#### Stage 1 gridsearch

```
snakemake --profile profile/slurm/ --snakefile stage_1
```

#### Stage 1 QC
##### Fit selection

After running stage 1 copy number fits will have been generated for each sample. In most cases, multiple viable fits will have been selected for each sample and a semi-qualitative quality control process needs to be applied to select the best fitting solution or exclude a sample should no fit be good enough.

Follow the guide on [fit selection](resources/quality_control_guide.md) and profile [Smoothing](#smoothing) to perform quality control and update the `{project}_fit_QC_predownsample.tsv` file. Once this step has been performed and a single acceptable fit has been selected for each sample, proceed to stage 2.

#### Stage 2

Provided quality control and fit selection was performed correctly, stage 2 of the pipeline can be performed which refits all copy number profiles using downsampled BAM files and the selected fits from stage 1.

As before, confirm the pipeline is configured correctly by running with the `dry-run` mode.

```
snakemake -n --profile profile/slurm/ --snakefile stage_2
```
and if the previous step ran without error then run the following:
```
snakemake --profile profile/slurm/ --snakefile stage_2
```
It is good practise to perform the QC performed in stage 1 on the results of stage 2.

</details>

## Output files

The final output is three types files located in the results directory `../sWGS_fitting/{project}_{bin}kb/absolute_POST_down_sampling/abs_cn_rds/`
- `*_ds_abs_fits.tsv` - Tab-seperated file containing absolute copy number profile metadata and fitting information
- `*_ds_absCopyNumber.rds` - QDNASeqmod class object containing binned copy number data in `rds` format
- `*_ds_absCopyNumber_segTable.tsv` - A multisample segment table consisting of copy number segment coordinates and associated segment values as unrounded absolute copy number

## Further details

### Smoothing

A subset of samples may require smoothing of segments in order to be viable. Read the [guide](resources/smoothing_guide.md) provided here to select and update which samples require smoothing. Once the `samplesheet.tsv` has been updated with the new `smooth` values, rerun stage 1 or auto as appropriate.

```
snakemake --profile profile/slurm/ -F all --snakefile stage_1 or auto
```

### Downsampling only

Where ploidy and purity are provided from alternative sources, users may not need perform a gridsearch. Provided all samples have a valid `precPloidy` and `precPurity`, the `processPrecomputed.R` script can generate all the required output files to skip directly to read depth normalisation (via read downsampling), at which point Stage 2 or auto can be executed as previously described.

```
Rscript scripts/processPrecomputed.R
```
or for container users
```
$SIGCMD Rscript scripts/processPrecomputed.R
```

### Reference genome

This pipeline currently supports both hg19 and hg38 reference genomes by specifying the correct build in the configuration file `config/config.yaml`. Currently there are no checks in place for using files aligned to an unsupported reference genome or to check that reference genomes are correct.

### Aligned read format

The pipeline accepts either BAM or CRAM files as the initial input, though it currently does not support the use of both concurrently. Users should specify the filetype in the `config/config.yaml` as either BAM or CRAM and provide a reference genome matching the one used in CRAM generation using the reference parameter in the `config/config.yaml`. If using the containerised implementation, it is recommended to place the reference genome within the input file directory to reduce unnecessary directory binding.

Using the CRAM implementation currently involves the decompresssion of CRAM files to BAM format as the underlying Rsamtools functions in QDNAseq will not load from CRAM. As such, the CRAM implementation will generate a larger diskspace footprint than the BAM implementation. The pipeline will currently remove decompressed BAMs once their required outputs are generated. 

### Profile troubleshooting
For most users, the default parameters should work well but in certain instances, these values should be modifed. 

For example, if users are not generating a sufficent number of high quality absolute copy number profiles, setting `filter_underpowered` to "FALSE" will show a larger range of fits which, although statistically underpowered, could be correct for given sample.

Another use case would be samples where the general purity range is known, for example cell line or LCM data, where the expected purity is known to be ~1.0. As such setting `purity_min` to 0.95 would restrict ploidy/purity combinations to enforce this expectation.

## Authors

* Philip Smith (@phil9s)

## Citation

Please cite this pipeline using the publication ["The copy-number landscape of recurrent ovarian high grade serous carcinoma"](https://doi.org/10.1038/s41467-023-39867-7) Smith & Bradley et al. 2023 - _Nature Communications_ and/or the version controlled zenodo repository [10.5281/zenodo.10040893](https://zenodo.org/doi/10.5281/zenodo.10040893).

