#!/bin/bash

echo "Basecalling sample ${snakemake_wildcards[sample_name]}"

dorado \
basecaller ${snakemake_params[model]} ${snakemake_input[pod5_dir]} \
    --min-qscore ${snakemake_params[min_q]} \
    > ${snakemake_output[basecall_bam]} 2> ${snakemake_log[basecall_log]}

samtools fastq \
    -T "*" \
    -0 ${snakemake_output[basecall_fastq]} \
    ${snakemake_output[basecall_bam]}
