#!/bin/bash

# create a bed file for each nucleotide
modkit motif bed \
    ${snakemake_input[ref]} \
    ${snakemake_wildcards[nt]} 0 > ${snakemake_output[temp_bed]}

# awk remove -ve strand
# sed add a column to indicate the nucleotide 
# awk replace column 4 with -

awk \
    '$6 != "-"' ${snakemake_output[temp_bed]} | \
sed \
    -e "s/$/\t${snakemake_wildcards[nt]}/" | \
awk \
    '{FS=OFS="\t"}{$4 = "-"}1' > \
    ${snakemake_output[modkit_out]}
