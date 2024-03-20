# Snakemake workflow: `SNP_GATK`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)



A Snakemake workflow for `calling snp`


## Getting started

### Install Snakemake

```bash
mamba create -c conda-forge -c bioconda -n env_name snakemake
mamba activate env_name
pip install snakemake-executor-plugin-slurm
```

### Get the workflow

```bash
git clone https://github.com/r1cheu/snp_gatk.git
cd snp_gatk
sudo singularity build snp_gatk.sif snp_call.def
### we alse provide a pre-built singularity image
```

