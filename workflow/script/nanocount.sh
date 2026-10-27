#!/bin/bash

echo "nanocount sample ${snakemake_wildcards[sample_name]}"

NanoCount \
    -i ${snakemake_input[aln_sorted_bam]} \
    -o ${snakemake_output[temp_out]}

tail -n+2 ${snakemake_output[temp_out]} | \
awk -v sample="${snakemake_wildcards[sample_name]}" \
    'BEGIN{FS=OFS="\t"} {$1 = sample "_" $1}1' > ${snakemake_output[nanocount_tsv]}