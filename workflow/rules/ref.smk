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


rule make_intervals:
    """Make intervals for variant calling"""
    output:
        bed="<ref>/intervals/{ivl}.bed",
    params:
        ivl=lambda wildcards: "\n".join([repr(ivl) for ivl in intervals[wildcards.ivl]]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/make_intervals/<ref>/intervals/{ivl}.bed.benchmark.txt"
    log:
        "<logs>/make_intervals/<ref>/intervals/{ivl}.bed.log",
    threads: 1
    shell:
        """
        echo -e "{params.ivl}" > {output.bed}
        """


rule make_all_intervals:
    """Make all intervals file for variant calling"""
    output:
        "<ref>/intervals.bed",
    params:
        ivl=lambda wildcards: "\n".join(
            [repr(x) for ivl in intervals.values() for x in ivl]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/make_all_intervals/<ref>/intervals.bed.benchmark.txt"
    log:
        "<logs>/make_all_intervals/<ref>/intervals.bed.log",
    threads: 1
    shell:
        """
        echo -e "{params.ivl}" > {output}
        """
