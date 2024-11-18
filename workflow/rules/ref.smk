rule genome_faidx:
    input:
        get_ref_path()
    output:
        get_ref_path() + ".fai"
    log:
        "logs/samtools/create_faidx.log"
    cache: True
    shell:
        "samtools faidx {input.ref} > {output} {log}"

rule genome_dict:
    input:
        get_ref_path()
    output:
        get_ref_name() + ".dict"
    log:
        "logs/samtools/create_dict.log"
    cache: True
    shell:
        "samtools dict {input} -o {output} > {log} 2>&1"

rule bwa_index:
    input:
        get_ref_path()
    output:
        multiext(get_ref_path(), ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    resources:
        mem_mb=4000
    cache: True
    shell:
        "bwa index {input} > {log} 2>&1"
