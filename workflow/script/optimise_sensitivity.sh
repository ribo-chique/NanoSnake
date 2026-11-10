#!/bin/bash

Rscript ${snakemake[scriptdir]}/optimise_sensitivity_from_calls.R \
    -i ${snakemake_input[calls_tsv]} \
    -o ${snakemake_output[accuracy_tsv]} \
    -a ${snakemake_params[optimise_accuracy]}