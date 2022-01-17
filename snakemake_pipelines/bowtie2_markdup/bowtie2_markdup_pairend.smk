(samples,reads,) = glob_wildcards("fastq/{sample}_{read}.fastq.gz")
reads=["R1","R2"]
samples=list(set(samples))
samtools_path="samtools/samtools-1.13/samtools"

rule all:
    input:
        expand('marked/{sample}.bam', sample=samples)


rule bowtie2:
    input:
        r1 = "fastq/{sample}_R1.fastq.gz",
    	r2 = "fastq/{sample}_R2.fastq.gz"
    output:
        temp("mapped/{sample}.bam")
    log:
        "logs/bowtie2/{sample}.log"
    params:
        bowtie2_path="Bowtie/bowtie2-2.4.2-linux-x86_64/bowtie2", 
        bowtie2_index="Homo_sapiens/UCSC/hg38/Sequence/Bowtie2_2.3.5.1_Index/hg38",
    threads: 8
    shell:
        "({params.bowtie2_path} -x {params.bowtie2_index} -p {threads} -1 {input.r1} -2 {input.r2}  | {samtools_path} view -Sb -@ {threads} > {output}) 2> {log}"

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
        "sambamba/sambamba-0.8.1-linux-amd64-static markdup -t {threads} {input} {output}"

