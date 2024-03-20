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
    shell:
        "samtools dict {input} > {output} > {log} 2>&1"

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
    shell:
        "bwa index {input} > {log} 2>&1"
