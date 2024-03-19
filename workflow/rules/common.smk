import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config/config.yaml"
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
units = pd.read_table(config["units"], dtype=str).set_index("sample", drop=False)


def is_single_end(sample):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[sample, "fq2"])


def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}.fastq.gz".format(**wildcards)


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[wildcards.sample, "platform"],
    )


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/dedup/{sample}.sorted.bam",
        sample=wildcards.sample,
    )
