samples, = glob_wildcards("bam/{sample}.bam")
samtools_path="/PATH/TO/SAMTOOLS/"

rule all:
    input:
        expand('marked/{sample}.bam', sample=samples)

rule index:
    input:
        "bam/{sample}.bam"
    output:
        temp("bam/{sample}.bam.bai")
    shell:
        "{samtools_path} index {input}"

rule sambamba_markdup:
    input:
        bai="bam/{sample}.bam.bai",
        bam="bam/{sample}.bam"
    output:
        "marked/{sample}.bam"
    threads: 4
    shell:
        "/PATH/TO/sambamba/sambamba-0.7.0-linux-static markdup -t {threads} {input.bam} {output}"

