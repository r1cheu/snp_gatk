#!/bin/bash
#SBATCH -o gatk.out
#SBATCH -p normal
#SBATCH -J Bench_gatk
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --mem 2000
snakemake --profile cluster --cache --rerun-incomplete
