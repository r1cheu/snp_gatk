rule call_variants:
    input:
        bam=get_sample_bams,
        bam_idx=get_sample_bams_idx,
        ref=get_ref_path(),
        fai=get_ref_path() + ".fai",
        idx=get_ref_name() + ".dict",
    output:
        gvcf=protected("results/called/gvcfs/{sample}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    resources:
        cpus_per_task=4,
        mem_mb=16000
    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' HaplotypeCaller "
        "-R {input.ref} -I {input.bam} -O {output.gvcf} "
        "-ERC GVCF > {log} 2>&1"

rule generate_sample_map:
    input:
        gvcfs=expand("results/called/gvcfs/{sample}.g.vcf.gz", sample=samples.index),
    output:
        temp("results/called/genome_db/sample_map")
    run:
        import os.path as osp
        with open(output[0], "w") as f:
            for gvcf in input.gvcfs:
                sample = osp.basename(gvcf).split(".")[0]
                f.write(f"{sample}\t{gvcf}\n")


rule genomics_db_import:
    input:
        "results/called/genome_db/sample_map",
    output:
        db=directory("results/called/genome_db/{chrom}_db"),
    log:
        "logs/gatk/genomics_db/genomics_db_import_{chrom}.log",
    benchmark:
        "results/called/genome_db/genomics_db_import_{chrom}.benchmark",
    resources:
        cpus_per_task=4,
        mem_mb=8000,
        slurm_partition='fat',
    shell:
        """
        gatk --java-options '-Xmx{resources.mem_mb}m' GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --batch-size 100 --reader-threads {resources.cpus_per_task} \
        --sample-name-map {input} \
        -L {wildcards.chrom} > {log} 2>&1
        """


rule joint_calling:
    input:
        db="results/called/genome_db/{chrom}_db",
        ref=get_ref_path(),
        idx=get_ref_name() + ".dict",
    output:
        vcf=protected("results/vcf/all.{chrom}.vcf.gz"),
    log:
        "logs/gatk/joint_calling/{chrom}_joint_calling.log",
    benchmark:
        "logs/gatk/joint_calling/{chrom}_joint_calling.benchmark",
    resources:
        mem_mb=32000,
        slurm_partition='fat',

    shell:
        "gatk --java-options '-Xmx{resources.mem_mb}m' GenotypeGVCFs "
        "-R {input.ref} -V gendb://{input.db} "
        "-O {output.vcf} > {log} 2>&1"

