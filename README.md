# sWGS-absoluteCN pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.10.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Philip Smith (@phil9s)

## Usage

Generate absolute copy number profiles from shallow whole genome sequencing data using a read depth normalised and allele frequency-anchored approach.

### Step 1: Clone the repo

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system.

```
wget https://github.com/Phil9S/sWGS-absoluteCN.git
cd sWGS-absoluteCN/
```

### Step 2: Install conda

Run the following to install conda whilst following the on-screen instructions.
- When asked to run `conda init` and initialise conda please respond with 'yes'

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $HOME/miniconda/
rm Miniconda3-latest-Linux-x86_64.sh
```

See installing [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for more information.

#### For those with Conda already installed

For systems where conda is already available the following requirements need to be met:
- conda must be available on the PATH
- conda version `4.8.3' or greater
- the location of the installation folder is required

Check the installed version of conda using the following:
```
conda -V
```
*If this command does not work then conda is also not available on the PATH*

Find your installation directory using the following:
```
whereis conda | sed 's%condabin/conda%%'
```

### Step 3: Installing additional dependencies

From within the repository directory, run the `install_env.sh` script to generate a conda environment and install custom packages:

```
install_env.sh $HOME/miniconda/'
```

If you used a previously installed conda build please use the conda or miniconda installation directory when running this section instead of '$HOME/miniconda/' to correctly initialise the conda environment.

*To be replaced with a built-in snakemake solution once possible*

### Step 4: Preparing the input files

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

#### slurm config.yaml

The slurm config.yaml (`profile/slurm/config.yaml`) contains the necessary information to configure the job submission parameters passed to `sbatch` on slurm-managed cluster enviroments. This includes the number of cocurrently sumbitted jobs, slurm account name, slurm partition name, and default job resources (though these are low and should work on almost any cluster).

### Updating the pipeline configuration

Environment-specific and pipeline-specific parameters need to be set for each run of this pipeline. While it is possible to manually edit the YAML files, a script has been provided to update the most frequently altered parameters programmatically.

With the `swgs-abscn` conda environment active, run the following code:

```
python update_configs.py
```

### Step 5: Stage_1

Once the pipeline and cluster parameters have been set and the samplesheet is prepared, the first stage of the pipeline is ready to run.

With the `swgs-abscn` conda environment active, run the following:

Confirm the pipeline is configured correctly, run the first stage using the `dry-run` mode.

```
snakemake -n --profile config/slurm/ --snakefile stage_1
```

If the previous step ran without error then run the following:

```
snakemake --profile config/slurm/ --snakefile stage_1
```

### Step 6: Stage 1 - QC1 and fit selection

### Step 7: Stage 2 - Downsampling and QC2

### Step 8: Stage 3 - Cohort-level filtering

## Quality control methodologies

### Quality control (smoothing)

### Quality control (fit selection)

