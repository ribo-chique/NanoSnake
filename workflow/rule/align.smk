rule bwa_index:
    input:
        ref=config["reference"]
    output:
        multiext(config["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    params:
        idx=config["bwa"]["idx-param"]
    conda:
        "../env/bwa.yaml"
    script:
        "../script/bwa_index.sh"


rule bwa_mem:
    input:
        basecall_fastq="results/01_basecall/{sample_name}/{sample_name}.basecall.fastq.gz",
        ref=config["reference"],
        index=rules.bwa_index.output
    output:
        multiext("results/02_align/{sample_name}/{sample_name}", aln_sorted_bam=".aln.sorted.bam", aln_fastq=".aln.fastq.gz")
    params:
        mem=config["bwa"]["mem-param"]
    conda:
        "../env/bwa.yaml"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=2,
        mem_mb_per_cpu=8000
    script:
        "../script/align.sh"
