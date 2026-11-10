def override_threshold_param(wildcards):
    override_threshold = config["pileup"]["override_threshold"]
    # if field is non-empty and it's not "default", use entire field
    if override_threshold and override_threshold != "default":
        per_base = override_threshold.split()
        return " ".join(f"--filter-threshold {base}" for base in per_base)
    # if field is empty or if it's "default", return nothing
    elif override_threshold == "default" or not override_threshold:
        return ""

def mod_threshold_param(wildcards):
    mod_threshold = config["pileup"]["mod_threshold"]
    # if field is non-empty, use entire field
    if mod_threshold:
        per_mod = mod_threshold.split()
        return " ".join(f"--mod-threshold {mod}" for mod in per_mod)
    # if field is empty, return nothing
    return ""

def pileup_region_params(wildcards):
    # if there are chrom specified, output region param
    if QC_CHROM != "all":
        return f"--region {wildcards.chrom}"
    # if field is empty, return nothing
    return ""

def pileup_input_files(wildcards):
    sample_bam = f"results/02_align/{wildcards.sample_name}/{wildcards.sample_name}.aln.sorted.bam"
    accuracy_tsv = ["results/04_sample_probs/{optimise_sample}/{optimise_sample}_optimise_accuracy.tsv".format(optimise_sample=optimise_sample) for optimise_sample in OPTIMISE_SAMPLE]

    # if there is no override param, then must provide threshold_tsv
    if not config["pileup"]["override_threshold"]:
        return {
            "aln_sorted_bam" : sample_bam,
            "accuracy_tsv": accuracy_tsv
        }
    
    # if there is override param, then don't need threshold_tsv
    else:
        return {
            "aln_sorted_bam" : sample_bam,
            "accuracy_tsv": ""
        }

if QC_CHROM != "all":
    rule pileup_chrom:
        input:
            unpack(pileup_input_files)
        output:
            pileup_bed = "results/05_pileup/{sample_name}/{chrom}/{sample_name}_{chrom}.pileup.bed"
        params:
            ref = config["reference"],
            region = pileup_region_params,
            other_params = config["pileup"]["pileup_params"],
            bed = f"results/fa2bed/{ref_name}.bed",
            include_bed = config["pileup"]["include_bed"],
            chrom_specific_bool = config["pileup"]["use_chrom_specific_threshold"],
            override_threshold = override_threshold_param,
            mod_threshold = mod_threshold_param,
            log_path = "results/05_pileup/{sample_name}/{chrom}/{sample_name}_{chrom}.pileup.log"
        resources:
            #24h
            runtime=1440,
            cpus_per_task=2,
            mem_mb=16000
        script:
            "../script/pileup.sh"

    rule concatenate_pileup:
        input:
            pileup_bed = expand("results/05_pileup/{sample_name}/{chrom}/{sample_name}_{chrom}.pileup.bed", chrom = QC_CHROM, allow_missing=True)
        output:
            pileup_bed = "results/05_pileup/{sample_name}/{sample_name}.pileup.bed"
        shell:
            "cat {input.pileup_bed} > {output.pileup_bed}"

else:
    rule pileup_all:
        input:
            unpack(pileup_input_files)
        output:
            pileup_bed = "results/05_pileup/{sample_name}/{sample_name}.pileup.bed"
        params:
            ref = config["reference"],
            region = pileup_region_params,
            other_params = config["pileup"]["pileup_params"],
            bed = f"results/fa2bed/{ref_name}.bed",
            include_bed = config["pileup"]["include_bed"],
            chrom_specific_bool = config["pileup"]["use_chrom_specific_threshold"],
            override_threshold = override_threshold_param,
            mod_threshold = mod_threshold_param,
            log_path = "results/05_pileup/{sample_name}/{sample_name}.pileup.log"
        resources:
            #24h
            runtime=1440,
            cpus_per_task=2,
            mem_mb=16000
        script:
            "../script/pileup.sh"

