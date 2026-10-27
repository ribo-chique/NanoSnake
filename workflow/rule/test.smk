rule basecall:

    input:
        config['input_dir']+'/{sample}'
    output:
        'result/01_basecall/{sample}/{sample}.basecall.fastq.gz'
    script:
        'workflow/script/basecall.sh'
