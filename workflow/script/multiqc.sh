#!/bin/bash

echo "Aggregating sample QC"

echo "Running bam2stats to all.stats"
# Run bam2stat then aggregate to AlignmentQC file - general stats in overall bam
python ${snakemake[scriptdir]}/bam2stats.py \
    ${snakemake_input[aln_sorted_bam]} > \
    ${snakemake_output[bam2stat_all]}

echo "Set up alnQC_mqc"
## Generate multiqc template
echo '
# id: alnQC
# plot_type: "table"
# section_name: "Alignment QC"
# file_format: "tsv"
 ' > ${snakemake_output[bam2stat_mqc]}

### Read one overall bam2stat output
### Get first line
### remove hash character, append line to template
head -n1 ${snakemake_output[bam2stat_all]} | \
sed 's/#//g' >> ${snakemake_output[bam2stat_mqc]}

echo "Aggregate bam2stats to alnQC_mqc"
### Read bam2stat output
### Read all lines that do not have a hash character
### Assign tab separators, remove extension, remove parent directory name, append line to template
grep -v '#' ${snakemake_output[bam2stat_all]} | \
awk 'BEGIN{FS=OFS="\t"} {gsub("\\.aln\\.sorted\\.bam","",$1)}{gsub("(.*/)+","",$1)}1' >> ${snakemake_output[bam2stat_mqc]}

echo "Set up nanocount_mqc"
# Summarise counts in Nanocount files
## Generate multiqc template
echo '
# id: NanoCount
# plot_type: "table"
# section_name: "NanoCount assigned read counts"
# file_format: "tsv"
File name	Proportion	Estimated count	TPM
' > ${snakemake_output[nanocount_mqc]}

echo "Aggregate nanocount to nanocount_mqc"
### Read individual nanocount output, append line to template
cat ${snakemake_input[nanocount_tsv]} >> ${snakemake_output[nanocount_mqc]}

echo "Set up counts_mqc"
# Run count then aggregate to AlignmentQC file - general stats in overall bam
## Generate multiqc template
echo '
# id: Assigned reads
# plot_type: "table"
# section_name: "Assigned counts"
File name	Counts
' > ${snakemake_output[samtools_countstat_mqc]}

echo "Aggreagate count.stats to counts_mqc"
### Read individual samtools counts output
### append line to template
cat ${snakemake_input[count_stat]} >> ${snakemake_output[samtools_countstat_mqc]}

echo "Generate multiQC"
# Generate multiqc report
multiqc \
    ${snakemake_params[output_dir]} \
    -f \
    --cl-config \
'extra_fn_clean_exts:
- "." ' \
    -n \
    ${snakemake_output[mqc_final]}
