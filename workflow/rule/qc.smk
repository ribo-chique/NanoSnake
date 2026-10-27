rule split_bam:
    input:
        aln_sorted_bam = "results/02_align/{sample_name}/{sample_name}.aln.sorted.bam"
    output:
        temp(multiext("results/02_align/{sample_name}/{sample_name}_{chrom}", aln_sorted_bam=".aln.sorted.bam", aln_fastq=".aln.fastq.gz", bam_idx=".aln.sorted.bam.bai"))
    params:
        chrom=config["QC"]["chrom"]
    resources:
        #24h
        runtime=1440,
        mem_mb=16000
    conda:
        "../env/multiqc.yaml"
    script:
        "../script/bamsplit.sh"

rule nanocount:
    input:
        aln_sorted_bam = "results/02_align/{sample_name}/{sample_name}.aln.sorted.bam"
    output:
        temp_out = temp("results/03_QC/{sample_name}/{sample_name}.temp"),
        nanocount_tsv = temp("results/03_QC/{sample_name}/{sample_name}_all.nanocount.tsv")
    conda:
        "../env/multiqc.yaml"
    resources:
        cpus_per_task=1,
        mem_mb=32000
    script:
        "../script/nanocount.sh"

# Find input files for the main QC step. common.smk defined chrom to be "all" + chroms in config.
def QC_input_files(wildcards):
    file_name = f"results/02_align/{wildcards.sample_name}/{wildcards.sample_name}"
    if wildcards.chrom == "all":
        return {
            "aln_sorted_bam": f"{file_name}.aln.sorted.bam",
            "aln_fastq": f"{file_name}.aln.fastq.gz"
            }
    else:
        return {
            "aln_sorted_bam": f"{file_name}_{wildcards.chrom}.aln.sorted.bam",
            "aln_fastq": f"{file_name}_{wildcards.chrom}.aln.fastq.gz"
            }

rule QC:
    input:
        unpack(QC_input_files)
    output:
        temp(multiext("results/03_QC/{sample_name}/{sample_name}_{chrom}", nanoq=".nanoq", nanostat=".nanostat", assigned=".assigned", count_stat=".count.stats"))
    conda:
        "../env/multiqc.yaml"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=4,
        mem_mb=64000
    script:
        "../script/qc.sh"

rule multiQC:
    input:
        # despite these being temp files, I'm listing them here to ensure the previous step runs first
        # excluding IVT for nanocount because it takes too much mem
        # label sample_chrom_bam_idx as input just in case these are needed, also only label the chrom that are not "all"
        nanocount_tsv = expand("results/03_QC/{sample_name}/{sample_name}_all.nanocount.tsv", sample_name = BIO_SAMPLE),
        aln_sorted_bam = expand("results/02_align/{sample_name}/{sample_name}.aln.sorted.bam", sample_name = SAMPLE),
        chrom_bam_idx = expand("results/02_align/{sample_name}/{sample_name}_{chrom}.aln.sorted.bam.bai", sample_name = SAMPLE, chrom = config["QC"]["chrom"]),

        nanoq = expand("results/03_QC/{sample_name}/{sample_name}_{chrom}.nanoq", sample_name = SAMPLE, chrom = CHROM),
        nanostat = expand("results/03_QC/{sample_name}/{sample_name}_{chrom}.nanostat", sample_name = SAMPLE, chrom = CHROM),
        count_stat = expand("results/03_QC/{sample_name}/{sample_name}_{chrom}.count.stats", sample_name = SAMPLE, chrom = CHROM)
    output:
        bam2stat_all = temp("results/03_QC/all.stats"),
        bam2stat_mqc = "results/03_QC/alnQC_mqc.txt",
        nanocount_mqc = "results/03_QC/nanocount_mqc.txt",
        samtools_countstat_mqc = "results/03_QC/counts_mqc.txt",
        mqc_final = "results/03_QC/multiqc.html"
    params:
        output_dir = "results/03_QC"
    conda:
        "../env/multiqc.yaml"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=1,
        mem_mb=4000
    script:
        "../script/multiqc.sh"