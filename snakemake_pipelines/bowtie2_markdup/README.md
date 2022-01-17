# Alignment and Mark Duplicates

This [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline performs alignment and removal of PCR duplicates from a fastq file.

# Output

A sorted `.bam` file for each aligned sampled with duplicates marked.

# Installation

Make sure to install and properly run snakemake, clone this repository and run the snakefile.

# Usage

All `.fastq` files must be inside a folder called `fastq`. After that run:

```
snakemake --snakefile /PATH/TO/SNAKEFILE/bowtie_markdup.smk --cores 60 
```

For paired-end reads, run:
```
snakemake --snakefile /PATH/TO/SNAKEFILE/bowtie2_markdup_pairend.smk --cores 60 
```


# Dependencies

- snakemake
- sambamba

# Pipeline steps

This pipelines uses `bowtie2` to align reads to hg38 and `sambamba` to mark duplicates.

