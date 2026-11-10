#!/bin/bash

threshold_tsv=${snakemake_input[accuracy_tsv]}
chrom_specific=${snakemake_params[chrom_specific_bool]}
override_threshold=${snakemake_params[override_threshold]}

# -n = if non-empty, -z = if empty
# if threshold_tsv exist and override is empty
if [[ -n "$threshold_tsv" && -z $override_threshold ]]; then

    filter_threshold=()
    # then read threshold tsv
    while IFS=$'\t' read -r CHROM NT RAW_AC DEFAULT_PROB DEFAULT_AC FINAL_PROB FINAL_AC
    do
        # while reading, if we are not using chrom specific threshold, whole table is used
        if [[ $chrom_specific == "False" ]]; then
            filter_threshold+=("--filter-threshold" "${NT}:${FINAL_PROB}")
        
        # if we are using chrom specific threshold, select for rows where chrom column matches with current region param
        # -n "${snakemake_wildcards[chrom]:-}" make sure wildcard is present and non-empty, guard against if QC_chrom == all and chrom_specific == True, although threshold_tsv should be unset
        elif [[ $chrom_specific == "True" && ${CHROM} == ${snakemake_wildcards[chrom]} && -n "${snakemake_wildcards[chrom]:-}" ]]; then
            filter_threshold+=("--filter-threshold" "${NT}:${FINAL_PROB}")
        fi
    done < <(tail -n+2 $threshold_tsv)

# if threshold_tsv does not exist and override param is non-empty
elif [[ -z "$threshold_tsv" && -n $override_threshold ]]; then

    # then override must be custom thresholds, use verbatim
    filter_threshold=($override_threshold)

# if threshold_tsv or override param is empty
elif [[ -z "$threshold_tsv" || -z $override_threshold ]]; then

    # then use default threshold, remain empty
    filter_threshold=()
fi

include_bed_bool=${snakemake_params[include_bed]}
# if pileup filters for bed file
bed_param=""
if [[ include_bed_bool == "True" ]]; then

    # then use --include-bed paramter
    bed_param="--include-bed ${snakemake_params[bed]}"
fi

# ${filter_threshold[@]+"$filter_threshold[@]}"} = if filter_threshold is non-empty, expand all elements of the array (ie --filter-threshold <...> --filter-threshold <...> ),
# if filter_threshold is empty then output nothing, which is not the same as an empty array

modkit pileup \
    --reference ${snakemake_params[ref]} \
    ${snakemake_params[region]} \
    ${snakemake_params[other_params]} \
    ${bed_param} \
    ${filter_threshold[@]+"${filter_threshold[@]}"} \
    ${snakemake_params[mod_threshold]} \
    --threads ${snakemake_resources[cpus_per_task]} \
    --log-filepath ${snakemake_params[log_path]} \
    ${snakemake_input[aln_sorted_bam]} \
    ${snakemake_output[pileup_bed]}
