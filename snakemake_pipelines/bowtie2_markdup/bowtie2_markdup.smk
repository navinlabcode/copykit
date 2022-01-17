samples, = glob_wildcards("fastq/{sample}.fastq.gz")
samtools_path="samtools-1.13/samtools"

rule all:
    input:
        expand('marked/{sample}.bam', sample=samples)


rule bowtie2:
    input:
        r1 = "fastq/{sample}.fastq.gz",
    output:
        temp("mapped/{sample}.bam")
    log:
        "logs/bowtie2/{sample}.log"
    params:
        bowtie2_path="bowtie2-2.4.4/bowtie2", 
        bowtie2_index="Homo_sapiens/UCSC/hg38/Sequence/Bowtie2_2.4.4_Index/hg38",
    threads: 8
    shell:
        "({params.bowtie2_path} -x {params.bowtie2_index} -p {threads} -U {input.r1}  | {samtools_path} view -Sb -@ {threads} > {output}) 2> {log}"

rule sort:
    input: 
        "mapped/{sample}.bam"
    output:
        temp("sort/{sample}.bam")
    threads: 4
    shell:
        "{samtools_path} sort {input} -@ {threads} -o {output}"

rule index:
    input:
        "sort/{sample}.bam"
    output:
        temp("sort/{sample}.bam.bai")
    shell:
        "{samtools_path} index {input}"

rule sambamba_markdup:
    input:
        "sort/{sample}.bam",
    output:
        "marked/{sample}.bam"
    threads: 4
    shell:
        "sambamba/sambamba-0.7.0-linux-static markdup -t {threads} {input} {output}"

