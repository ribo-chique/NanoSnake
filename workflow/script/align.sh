#!/bin/bash

echo "Aligning sample ${snakemake_wildcards[sample_name]}"

bwa mem \
    ${snakemake_params[mem]} \
    -t ${snakemake_resources[cpus_per_task]} \
    ${snakemake_input[ref]} ${snakemake_input[basecall_fastq]} | \
samtools sort \
    -@ ${snakemake_resources[cpus_per_task]} \
    -o ${snakemake_output[aln_sorted_bam]} -


samtools index \
    ${snakemake_output[aln_sorted_bam]}

samtools fastq \
    -T "*" \
    -0 ${snakemake_output[aln_fastq]} \
    ${snakemake_output[aln_sorted_bam]}
