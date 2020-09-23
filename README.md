# Snakemake workflow: sWGS-absoluteCN

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Philip Smith (@phil9s)

## Usage

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach.

### Step 1: Install conda

Run the following to install conda

```
wget -c -q 
bash Miniconda3-latest-Linux-x86_64.sh
```
Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

Make a note of where the miniconda directory is installed as this directory location is required later.

### Step 2: Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
wget https://github.com/Phil9S/sWGS-absoluteCN.git
cd sWGS-absoluteCN/
```
### Step 3: Installing additional dependencies

From within the repository directory, run the `install_env.sh` script to generate a conda environment and install custom packages*.

```
install_env.sh /home/miniconda_directory/
```

Use the conda/miniconda directory noted earlier when running this section to correctly initialise the conda environment.

*To be replaced with a built-in snakemake solution once possible*

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

#### slurm config.yaml


#### Config 

### Step 5: Stage_1

With the `swgs-abscn` conda environment active, run the following code:

```
snakemake --profile config/slurm/ --snakefile stage_1
```

### Step 6: Stage 1 - QC1 and fit selection

### Step 7: Stage 2 - Downsampling and QC2

### Step 8: Stage 3 - Cohort-level filtering

## Quality control methodologies

### Quality control (smoothing)

### Quality control (fit selection)

