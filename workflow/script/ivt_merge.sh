#!/bin/bash

samtools merge \
    ${snakemake_input[aln_sorted_bam]} \
    -o ${snakemake_output[aln_sorted_bam]}

samtools index ${snakemake_output[aln_sorted_bam]}
