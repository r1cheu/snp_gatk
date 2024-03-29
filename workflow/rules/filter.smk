rule get_merge_list:
    input: 
        expand("results/vcf/all.{chrom}.vcf.gz", chrom=CHROMOSOMES),
    output:
        "results/vcf/merge_list.txt",
    run:
        with open(output[0], "w") as f:
            for inp in input:
                f.write(inp + '\n')


rule merge_vcfs:
    input: 
        "results/vcf/merge_list.txt"
    output:
        "results/all.vcf.gz"
    log:
        "logs/gatk/merge_vcfs.log"
    benchmark:
        "logs/gatk/merge_vcfs.benchmark"
    shell: 
        "gatk MergeVcfs -I {input} -O {output} > {log} 2>&1"

rule select_snp:
    input:
        "results/all.vcf.gz"
    output:
        "results/all.snp.vcf.gz"
    log:
        "logs/gatk/select_snp.log"
    benchmark:
        "logs/gatk/select_snp.benchmark"
    shell:
        "gatk SelectVariants -V {input} -select-type SNP -O {output} > {log} 2>&1"


rule filter_snp:
    input:
        vcf="results/all.snp.vcf.gz",
        ref="resources/IRGSP-1.0_genome.fasta"
    output:
        "results/all.snp.filter.vcf.gz"
    log:
        "logs/gatk/filter_snp.log"
    benchmark:
        "logs/gatk/filter_snp.benchmark"
    shell:
        "gatk VariantFiltration -R {input.ref} --variant {input.vcf} "
        "--cluster-size 3 --cluster-window-size 10 --filter-expression "
        "'QD<10.00' --filter-name lowQD --filter-expression 'FS>15.000' "
        "--filter-name highFS --genotype-filter-expression 'DP>200||DP<5' "
        "--genotype-filter-name InvalidDP --output {output} > {log} 2>&1"
    

rule final_vcf:
    input:
        "results/all.snp.filter.vcf.gz"
    output:
        "results/all.snp.final.vcf.gz"
    log:
        "logs/gatk/final_vcf.log"
    shell:
        "bcftools view -f PASS {input} -Oz -o {output}"
