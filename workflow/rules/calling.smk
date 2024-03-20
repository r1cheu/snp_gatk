rule call_variants:
    input: 
        bam=get_sample_bams,
        ref="resources/IRGSP-1.0_genome.fasta",
        idx="resources/IRGSP-1.0_genome.dict",
    output: 
        gvcf=protected("results/called/{sample}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    run:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.gvcf} "
        "-ERC GVCF"


rule genomics_db_import:
    input: 
        gvcfs=expand("results/called/{sample}.g.vcf.gz", sample=samples.index),
    output: 
        db="results/called/genomics_db",
    log:
        "logs/gatk/genomics_db_import.log",
    threads: 5,
    run:
        "gatk GenomicsDBImport --genomicsdb-workspace-path {output.db} "
        "--batch-size 100 --reader-threads {threads} {input.gvcfs}"


rule joint_calling:
    input: 
        db="results/called/genomics_db",
        ref="resources/IRGSP-1.0_genome.fasta",
        idx="resources/IRGSP-1.0_genome.dict",
    output: 
        vcf="results/called/joint.vcf.gz",
    log:
        "logs/gatk/joint_calling.log",
    run:
        "gatk GenotypeGVCFs -R {input.ref} -V gendb://{input.db} "
        "-O {output.vcf} 2> {log}"
