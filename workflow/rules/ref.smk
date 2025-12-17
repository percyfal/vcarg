rule picard_create_sequence_dictionary:
    """Create sequence dictionary"""
    output:
        "{prefix}.dict",
    input:
        "{prefix}.fasta",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/{prefix}.dict.benchmark.txt"
    log:
        "logs/{prefix}.dict.log",
    threads: 1
    shell:
        """
        picard CreateSequenceDictionary --REFERENCE {input} > {log} 2>&1
        """
