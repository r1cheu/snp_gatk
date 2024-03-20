import os.path as osp

rule fastp_trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("results/trimmed/{sample}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}.2.fastq.gz")
    log:
        html="logs/fastp/{sample}.html", json="logs/fastp/{sample}.json"
    resources:
        ntasks=4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {log.html} -j {log.json} -w {resources.ntasks}"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.sorted.bam"),
    log:
        "logs/bwa/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort_order="coordinate",
    resources:
        ntasks=4
    shell:
        "bwa mem -t {resources.ntasks} {params.extra} {params.index} {input.reads} | "
        "samtools sort -T {output}.tmp -o {output} - > {log} 2>&1"
    

rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics=temp("results/qc/dedup/{sample}.metrics")
    log:
        "logs/picard/dedup/{sample}.log" 
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} "
        "--VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 "
        "--ASSUME_SORT_ORDER 'coordinate' --CREATE_INDEX true > {log} 2>&1"
    