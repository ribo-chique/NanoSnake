#!/bin/bash

Rscript ${snakemake[scriptdir]}/optimise_sensitivity_from_calls.R \
    -i ${snakemake_input[out_dir]} \
    -a ${snakemake_params[optimise_accuracy]} \
    -p ${snakemake_resources[cpus_per_task]}