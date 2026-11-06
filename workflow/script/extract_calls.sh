#!/bin/bash

# sample mod calls for each input and chrom
modkit extract calls \
    --reference ${snakemake_params[ref]} \
    --preload-references \
    ${snakemake_params[other_params]} \
    --threads ${snakemake_resources[cpus_per_task]} \
    ${snakemake_params[output_params]} \
    ${snakemake_input[aln_sorted_bam]} \
    ${snakemake_output[calls_tsv]}