rule genome_faidx:
    input: 
        "resources/IRGSP-1.0_genome.fasta"
    output:
        "resources/IRGSP-1.0_genome.fasta.fai"
    cache: True
    shell:
        "samtools faidx {input}"

rule genome_dict:
    input: 
        "resources/IRGSP-1.0_genome.fasta"
    output:
        "resources/IRGSP-1.0_genome.dict"
    log:
        "logs/samtools/create_dict.log"
    cache: True
    run:
        "samtools dict {input} > {output} 2> {log}"

rule bwa_index:
    input: 
        "resources/IRGSP-1.0_genome.fasta" 
    output: 
        multiext("resources/IRGSP-1.0_genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    resources:
        mem_mb=3000
    cache: True
    run:
        "bwa index {input} 2> {log}"
