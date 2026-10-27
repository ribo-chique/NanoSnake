rule tar_unpack:

    input:
        tgz=config["input_dir"]/{sample_name}_blow5.tgz"
    output:
        blow5_dir="config["input_dir"]/{sample_name}_blow5"
    resources:
        # 24h
        runtime: 1440 
        mem_mb: 4000
    shell:
        """
        tar -xvzf {input.tgz}
        """

rule blow5_to_pod5:

    input:
        blow5_dir="config["input_dir"]/{sample_name}_blow5/"
    output:
        pod5="result/01_pod5/{sample_name}/{sample_name}.pod5"
    log:
        "result/01_pod5/{sample_name}/{sample_name}.log"
    resources:
        # 24h
        runtime: 1440 
        threads: 16
        mem_mb: 16000
    shell:
        """
        blue-crab s2p {input.blow5_dir} -o {output.pod5} -p {resourse.threads} -v
        """
