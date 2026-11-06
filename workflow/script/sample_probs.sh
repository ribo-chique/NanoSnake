#!/bin/bash

# sample mod calls for each input and chrom
modkit sample-probs \
    --reference ${snakemake_params[ref]} \
    --preload-references \
    --include-bed ${snakemake_params[bed]} \
    --threads ${snakemake_resources[cpus_per_task]} \
    ${snakemake_params[config_params]} \
    ${snakemake_params[output_params]} \
    ${snakemake_input[aln_sorted_bam]}