#!/bin/bash

# concatenate input
cat ${snakemake_input[modkit_out]} > ${snakemake_output[bed]}

# sort text in column 1, then sort numbers in column 2
# -s stable sort, keep original order if column values are identical 
sort -sk1,1 -k2,2n -o ${snakemake_output[bed]} ${snakemake_output[bed]}