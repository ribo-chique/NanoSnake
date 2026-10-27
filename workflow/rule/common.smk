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

print("samples:", SAMPLE)
print("input directory:", config["input_dir"])

CHROM = ["all"] + config["QC"]["chrom"]

STEPS = ["01_basecall","02_align","03_QC"]
EXTENSION = ["basecall.bam", "aln.sorted.bam"]

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

    return final_output
