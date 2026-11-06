def parse_thresholds(wildcards):
    if config["pileup"]["override_threshold"] == "default":
        return ("")
    elif config["pileup"]["override_threshold"] != []:

        return (f"--filter-threshold config["pileup"]["override_threshold"]")
    elif config["pileup"]["override_threshold"] == []:
        if config["pileup"]["use_chrom_specific_threshold"] == False
        elif config["pileup"]["use_chrom_specific_threshold"] == True

def parse_thresholds_inputs(wildcards):
    accuracy_tsv = "results/04_sample_probs/{sample_name}/{chrom}/optimise_accuracy.tsv"

rule parse_thresholds:
    input:


# Get name of reference file without extension
ref_name = os.path.splitext(os.path.basename(config["reference"]))[0]

# Convert reference fa to bed, for mod sampling 
rule fa2bed_init:
    input:
        ref = config["reference"]
    output:
        temp_bed = temp(f"results/fa2bed/{ref_name}_{{nt}}.bed2"),
        modkit_out = temp(f"results/fa2bed/{ref_name}_{{nt}}.bed")
    script:
        "../script/fa2bed_init.sh"

rule fa2bed_end:
    input:
        modkit_out = expand(f"results/fa2bed/{ref_name}_{{nt}}.bed", nt = ["A", "C", "G", "T"])
    output:
        bed = f"results/fa2bed/{ref_name}.bed"
    script:
        "../script/fa2bed_end.sh"

# Find input files for the sample-probs step. If ivt_merge was ran then add merged bam to input
def probs_input_files(wildcards):
    sample_bam = f"results/02_align/{wildcards.sample_name}/{wildcards.sample_name}.aln.sorted.bam"
    merge_bam = "results/02_align/IVT_merge/IVT_merge.aln.sorted.bam"

    if os.path.exists(merge_bam):
        return {
            "aln_sorted_bam": [sample_bam, merge_bam]
            }
    else:
        return {
            "aln_sorted_bam": sample_bam
            }

# Separate rules for when there are chroms specified or not, because snakemake can only take conditional function in input but not output nor params
def probs_output_params_chrom(wildcards):
    return [
            f"--out-dir results/04_sample_probs/{wildcards.sample_name}/{wildcards.chrom}",
            f"--region {wildcards.chrom}",
            f"--prefix {wildcards.sample_name}_{wildcards.chrom}",
            f"--log-filepath results/04_sample_probs/{wildcards.sample_name}/{wildcards.chrom}/{wildcards.sample_name}_{wildcards.chrom}.probs.log"
        ]

rule sample_probs_chrom:
    input:
        unpack(probs_input_files)
    params:
        ref = config["reference"],
        bed = f"results/fa2bed/{ref_name}.bed",
        config_params = config["sample_probs"]["params"],
        output_params = probs_output_params_chrom
    output:
        probs_tsv = "results/04_sample_probs/{sample_name}/{chrom}/{sample_name}_{chrom}_probabilities.tsv"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=4,
        mem_mb=16000
    script:
        "../script/sample_probs.sh"

def probs_output_params_all(wildcards):
        return [
            f"--out-dir results/04_sample_probs/{wildcards.sample_name}",
            f"--prefix {wildcards.sample_name}",
            f"--log-filepath results/04_sample_probs/{wildcards.sample_name}/{wildcards.sample_name}.probs.log"
        ]

rule sample_probs_all:
    input:
        unpack(probs_input_files)
    params:
        ref = config["reference"],
        bed = f"results/fa2bed/{ref_name}.bed",
        config_params = config["sample_probs"]["params"],
        output_params = probs_output_params_all
    output:
        probs_tsv = "results/04_sample_probs/{sample_name}/{sample_name}_probabilities.tsv"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=4,
        mem_mb=16000
    script:
        "../script/sample_probs.sh"

def extract_calls_sensitivity_params_chrom(wildcards):
    return [
            f"--region {wildcards.chrom}"
        ]

rule extract_calls_sensitivity_chrom:
    input:
        unpack(probs_input_files)
    params:
        ref = config["reference"],
        other_params = "--bgzf --num-reads 10000 --no-filtering --mapped-only",
        output_params = extract_calls_sensitivity_params_chrom,
    output:
        calls_tsv = "results/04_sample_probs/{sample_name}/{chrom}/{sample_name}_{chrom}_extract_calls.tsv.gz"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=2,
        mem_mb=16000
    script:
        "../script/extract_calls.sh"

rule sensitivity_chrom:
    input:
        calls_tsv = "results/04_sample_probs/{sample_name}/{chrom}/{sample_name}_{chrom}_extract_calls.tsv.gz",
        out_dir = "results/04_sample_probs/{sample_name}/{chrom}/"
    params:
        optimise_accuracy = config["optimise_sensitivity"]["desired_sensitivity"]
    output:
        accuracy_tsv = "results/04_sample_probs/{sample_name}/{chrom}/optimise_accuracy.tsv"
    conda:
        "../env/r-optimise.yaml"
    resources:
        #24h
        runtime=1440,
        cpus_per_task=2,
        mem_mb=16000
    script:
        "../script/optimise_sensitivity.sh"