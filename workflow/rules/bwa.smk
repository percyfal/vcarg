BWA_INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]


rule bwa_index:
    """Make bwa index"""
    output:
        index=expand("bwa_index/{{prefix}}.fasta.{sfx}", sfx=BWA_INDEX_SUFFIX),
    input:
        fasta="ref/{prefix}.fasta",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "benchmarks/bwa_index/{prefix}.fasta.benchmark.txt"
    log:
        "logs/bwa_index/{prefix}.fasta.log",
    threads: 1
    shell:
        """
        bwa index -p bwa_index/{wildcards.prefix}.fasta {input.fasta} > {log} 2>&1
        """


rule bam_bwa_srr:
    """Map SRR data on the fly and convert to bam, without marking duplicates.

    NB: This will use the SampleName column to identify the sample
    (not the SRS id)."""
    output:
        bam=temp("bam_bwa_srr/{samplealias}/{srrun}.sort.bam"),
        bai=temp("bam_bwa_srr/{samplealias}/{srrun}.sort.bam.bai"),
    input:
        srr="sra/{srrun}/{srrun}.sra",
        reference=f"ref/{config['reference']}",
        index=expand(
            "{ref}.{sfx}", ref=f"bwa_index/{config['reference']}", sfx=BWA_INDEX_SUFFIX
        ),
    params:
        rg=get_read_group,
        options="-M",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "benchmarks/bam_bwa_srr/{samplealias}/{srrun}.sort.bam.benchmark.txt"
    log:
        "logs/bam_bwa_srr/{samplealias}/{srrun}.sort.bam.log",
    threads: 12
    priority: 50
    shell:
        """
        fasterq-dump -e {threads} --skip-technical -Z --split-spot {input.srr} |
            bwa mem {params.rg} -p -t {threads} {params.options} {input.reference} - 2>> {log}|
            samtools fixmate -m - /dev/stdout |
            samtools sort - |
            samtools view -h -b -o {output.bam} > {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """
