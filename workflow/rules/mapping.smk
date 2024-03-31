import os.path as osp


rule fastp_trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("results/trimmed/{sample}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}.2.fastq.gz")
    log:
        html="logs/fastp/{sample}.html", json="logs/fastp/{sample}.json"
    benchmark:
        "logs/fastp/{sample}.bench"
    resources:
        cpus_per_task=4,
        mem_mb=4000
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {log.html} -j {log.json} -w {resources.cpus_per_task}"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        "logs/bwa/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort_order="coordinate",
    resources:
        mem_mb=2000,
        cpus_per_task=4
    shell:
        "bwa mem -t {resources.cpus_per_task} -M {params.extra} {params.index} {input.reads} | "
        "samtools view -F 0x100 -Sb -o {output} - > {log} 2>&1"


rule sort_bam:
    input:
        "results/mapped/{sample}.bam"
    output:
        temp("results/mapped/{sample}.sorted.bam")
    log:
        "logs/picard/sort/{sample}.log"
    resources:
        mem_mb=8000
    shell:
        "gatk SortSam -I {input} -O {output} --SORT_ORDER coordinate > {log} 2>&1"


# 
rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics=temp("results/qc/dedup/{sample}.metrics")
    log:
        "logs/picard/dedup/{sample}.log" 
    benchmark:
        "logs/picard/dedup/{sample}.bench"
    resources:
        mem_mb=24000
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} "
        "--VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 "
        "--ASSUME_SORT_ORDER 'coordinate' > {log} 2>&1"


rule index_bam:
    input:
        "results/dedup/{sample}.sorted.bam"
    output:
        temp("results/dedup/{sample}.sorted.bam.bai")
    shell:
        "samtools index {input} {output}" 