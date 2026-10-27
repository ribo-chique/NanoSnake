rule bwa_map:
    input:
        pod5_dir="00_pod5/{sample_name}"
    output:
        multiext("result/01_basecall/{sample_name}/{sample_name}", basecall_bam=".basecall.bam", basecall_log=".basecall.log", basecall_fastq=".basecall.fastq")
    params:
        model=config["dorado"]["basecalling_model"],
        min_q=config["dorado"]["min_q"]
    log:
        "result/01_basecall/{sample_name}/{sample_name}.log"
    container:
        config["apptainer_path"]
    conda:
        "workflow/env/samtools.yaml"
    resources:
        slurm_partition="gpu",
        slurm_qos="gpu",
        #24h
        runtime=1440,
        gres="gpu:1",
        cpus_per_task=8,
        mem_mb_per_cpi=8000
    script:
        "script/basecall.sh"
