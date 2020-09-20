# Snakemake workflow: sWGS-absoluteCN

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Philip Smith (@phil9s)

## Usage

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach.

### Step 1: Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

### Step 2: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Installing additional dependencies

Run the following lines to add R packages which are unavailble from conda or are customised for this workflow to the conda environment path.

```
install_dependencies.sh
```

### Step 4: Preparing the input files

#### Sample sheet

The workflow requires a single input file `sample_sheet.tsv` which is a tab-separated document detailing sample names, patient names, smoothing booleans, TP53 allele frequencies, and BAM file locations. The table below demonstrates the basic schema required and the workflow will validate this file prior to running.

|PATIENT_ID|SAMPLE_ID|TP53freq|smooth|file         |
|----------|---------|--------|------|-------------|
|PAT1      |SAM1     |0.45    |TRUE  |data/SAM1.bam|
|PAT1      |SAM2     |0.55    |FALSE |data/SAM2.bam|

An example sample_sheet.tsv is included in this repository.

#### config.yaml

The config.yaml should be edited to contain the necessary information for the pipeline you wish to run. This includes the location of the samplesheet.tsv, bin size, project name, and output directory.

### Step 5: Copy number profile generate

With the conda environment run the following code:

```
snakemake -n --snakefile relativeQC --use-conda
```

### Step 6: Stage 1 - QC1 and fit selection


### Step 7: Stage 2 - Downsampling and QC2


### Step 8: Stage 3 - Cohort-level filtering absolute fit generation

## Quality control methodologies

### Quality control (smoothing)

### Qaulity control (fit selection)


