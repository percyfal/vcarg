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


rule samtools_faidx:
    """Index fasta file using samtools.faidx"""
    output:
        "{prefix}.fasta.fai",
    input:
        "{prefix}.fasta",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "<benchmarks>/samtools_faidx/{prefix}.fasta.fai.benchmark.txt"
    log:
        "<logs>/samtools_faidx/{prefix}.fasta.fai.log",
    threads: 1
    shell:
        """
        samtools faidx {input} > {log} 2>&1
        """
