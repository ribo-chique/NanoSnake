#!/bin/bash

echo "QC sample ${snakemake_wildcards[sample_name]}"

# Nanoq read stats
nanoq \
    -i ${snakemake_input[aln_fastq]} \
    -s \
    -vvv \
    -r ${snakemake_output[nanoq]}

# Nanostat read stats, similar to nanoq
NanoStat \
    --bam ${snakemake_input[aln_sorted_bam]} \
    --threads ${snakemake_resources[cpus_per_task]} \
    --name ${snakemake_output[nanostat]}


# Count number of mapped reads in overall or per gene bam, similar to NanoCount
# no secondary alignment, then cut output read_id and the chrom mapped
samtools view \
    -F 256 \
    -@ ${snakemake_resources[cpus_per_task]} \
    ${snakemake_input[aln_sorted_bam]} | \
cut -f 1,3 > \
    ${snakemake_output[assigned]} 

# count lines in file = number of primary mapped reads
# this is in the format "[number][space][file_name.extension]"
# sed remove extension
# awk assigns "[number][space][file_name]" as "[field1][space][field2]"
# then print "[field1][tab][field2]"

wc -l ${snakemake_output[assigned]} |\
sed 's/.assigned//g' | \
awk '{gsub("(.*/)+","",$2)}{print $2"\t"$1}' > \
    ${snakemake_output[count_stat]}