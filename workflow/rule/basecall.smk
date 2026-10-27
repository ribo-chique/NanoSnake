rule basecall:
    input:
        pod5_dir="00_pod5/{sample_name}/"
    output:
        multiext("results/01_basecall/{sample_name}/{sample_name}", basecall_bam=".basecall.bam", basecall_fastq=".basecall.fastq.gz")
    params:
        model=config["dorado"]["basecalling_model"],
        min_q=config["dorado"]["min_q"]
    log:
        basecall_log="results/01_basecall/{sample_name}/{sample_name}.log"
    container:
        config["apptainer_path"]
    conda:
        "../env/samtools.yaml"
    resources:
        slurm_partition="gpu",
        slurm_qos="gpu",
        #48h
        runtime=2880,
        gres="gpu:h100nvl:1",
        cpus_per_gpu=12,
        mem_mb=96000
    script:
        "../script/basecall.sh"
