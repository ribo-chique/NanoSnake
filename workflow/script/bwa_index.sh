#!/bin/bash

echo "Indexing ${snakemake_input[ref]}"

bwa index ${snakemake_params[idx]} ${snakemake_input[ref]}


