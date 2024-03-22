rule call_variants:
    input: 
        bam=get_sample_bams,
        bam_idx=get_sample_bams_idx,
        ref="resources/IRGSP-1.0_genome.fasta",
        fai="resources/IRGSP-1.0_genome.fasta.fai",
        idx="resources/IRGSP-1.0_genome.dict",
    output: 
        gvcf=protected("results/called/gvcfs/{sample}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    benchmark:
        "logs/gatk/haplotypecaller/{sample}.benchmark",
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.gvcf} "
        "-ERC GVCF > {log} 2>&1"


rule genomics_db_import:
    input: 
        gvcfs=expand("results/called/gvcfs/{sample}.g.vcf.gz", sample=samples.index),
    output: 
        db=directory("results/called/genome_db/{chrom}_db"),
    log:
        "logs/gatk/genomics_db/genomics_db_import_{chrom}.log",
    benchmark:
        "logs/gatk/genomics_db/genomics_db_import_{chrom}.benchmark",
    resources:
        cpus_per_task=4,
    shell:
        """
        gatk GenomicsDBImport --genomicsdb-workspace-path {output.db} \
        --batch-size 100 --reader-threads {resources.cpus_per_task} \
        -V {input.gvcfs} \
        -L {wildcards.chrom} > {log} 2>&1
        """


rule joint_calling:
    input: 
        db="results/called/genome_db/{chrom}_db",
        ref="resources/IRGSP-1.0_genome.fasta",
        idx="resources/IRGSP-1.0_genome.dict",
    output: 
        vcf="results/vcf/all.{chrom}.vcf.gz",
    log:
        "logs/gatk/joint_calling/{chrom}_joint_calling.log",
    benchmark:
        "logs/gatk/joint_calling/{chrom}_joint_calling.benchmark",
    shell:
        "gatk GenotypeGVCFs -R {input.ref} -V gendb://{input.db} "
        "-O {output.vcf} 2> {log} 2>&1"

