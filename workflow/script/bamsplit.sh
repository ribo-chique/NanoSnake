#!/bin/bash

samtools view \
    ${snakemake_input[aln_sorted_bam]} \
    ${snakemake_wildcards[chrom]} \
    -o ${snakemake_output[aln_sorted_bam]}

samtools index ${snakemake_output[aln_sorted_bam]}

samtools fastq \
    ${snakemake_output[aln_sorted_bam]} \
    -0 ${snakemake_output[aln_fastq]}