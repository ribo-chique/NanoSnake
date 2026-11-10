import glob
import os
import re
from snakemake.io import expand, glob_wildcards

SAMPLE = glob.glob(config["input_dir"]+"/**")
SAMPLE = [os.path.basename(path) for path in SAMPLE]

ivt_pattern = re.compile(r"ivt", re.IGNORECASE)

IVT_SAMPLE = [
    name for name in SAMPLE
    if ivt_pattern.search(name)
]
BIO_SAMPLE = [
    name for name in SAMPLE
    if not ivt_pattern.search(name)
]

SAMPLE_plus_merge = list(SAMPLE)
if len(IVT_SAMPLE) > 1:
    SAMPLE_plus_merge.append("IVT_merge")
    OPTIMISE_SAMPLE = ["IVT_merge"]
else:
    OPTIMISE_SAMPLE = IVT_SAMPLE

CHROM = ["all"] + config["QC"]["chrom"]

if config["QC"]["chrom"] == []:
    QC_CHROM = "all"
else:
    QC_CHROM = config["QC"]["chrom"]

if config["optimise_sensitivity"]["chrom"] == []:
    OPTIMISE_CHROM = QC_CHROM
else:
    OPTIMISE_CHROM = config["optimise_sensitivity"]["chrom"]

def get_final_output():
    
    final_output = expand(
        "results/01_basecall/{sample}/{sample}.basecall.bam",
        sample = SAMPLE
    )

    final_output.extend(
        expand(
        "results/02_align/{sample}/{sample}.aln.sorted.bam",
        sample = SAMPLE
    ))

    final_output.append(
        "results/03_QC/multiqc.html"
    )

    final_output.append(
        f"results/fa2bed/{ref_name}.bed"
    )

    if QC_CHROM == "all":
        final_output.extend(
            expand(
                "results/04_sample_probs/{sample}/{sample}_probabilities.tsv",
                sample = SAMPLE_plus_merge
                )
        )
    else:
        final_output.extend(
            expand(
            "results/04_sample_probs/{sample}/{chrom}/{sample}_{chrom}_probabilities.tsv",
            sample = SAMPLE_plus_merge,
            chrom = QC_CHROM
            )
        )
    
    final_output.extend(
        expand(
            "results/04_sample_probs/{sample}/{sample}_optimise_accuracy.tsv",
            sample = OPTIMISE_SAMPLE
            )
        )

    final_output.extend(
        expand(
            "results/05_pileup/{sample}/{sample}.pileup.bed",
            sample = SAMPLE
        )
    )
    return final_output

print("samples:", SAMPLE)
print("input directory:", config["input_dir"])
print("reference file:", config["reference"])
print("chroms to QC:", QC_CHROM)
print("IVT samples:", IVT_SAMPLE)