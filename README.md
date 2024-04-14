# Snakemake workflow: `SNP_GATK`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)



A Snakemake workflow for `calling snp`


## Getting started

### Install miniforge
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

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

### Prepare Samples
Note the workflow only accept suffix *.1.fastq.gz and *.2.fastq.gz. You need to
rename the fastq.gz file before start.
I will add a option to use arbitrary suffix. 

modify config/samples.tsv and config/units.tsv.
or use the provided python script
```python
python workflow/script/gen_input_csv.py -I /dir/path/to/fastq.gz -O config
```

### Modify the config file
modify the config/config.yaml file to fit your needs.
e.g. sif: /path/to/snp_gatk.sif
### Run the workflow

```bash
snakemake --profile cluster --cache
```