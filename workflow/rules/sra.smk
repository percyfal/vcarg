rule sra_prefetch:
    """Prefetch sra record. Note that output name is determined by
    prefetch."""
    output:
        srrun=temp("<project>/sra/{srrun}/{srrun}.sra"),
    conda:
        "../envs/sratools.yaml"
    benchmark:
        "benchmarks/<project>/sra_prefetch/sra/{srrun}.benchmark.txt"
    log:
        "logs/<project>/sra_prefetch/sra/{srrun}.log",
    priority: 0
    threads: 1
    shell:
        """
        prefetch {wildcards.srrun} -o {output.srrun} -p > {log} 2>&1
        """


rule fastq_fasterq_dump:
    """Dump fastq from sra file using fasterq-dump."""
    output:
        fastq="<project>/fastq_fasterq_dump/{srrun}/{srrun}.fastq.gz",
    input:
        sra="<project>/sra/{srrun}/{srrun}.sra",
    conda:
        "../envs/sratools.yaml"
    benchmark:
        "benchmarks/<project>/fastq_fasterq_dump/{srrun}/{srrun}.fastq.gz.benchmark.txt"
    log:
        "logs/<project>/fastq_fasterq_dump/{srrun}/{srrun}.fastq.gz.log",
    threads: 12
    shell:
        """
        fasterq-dump -e {threads} --skip-technical -Z --split-spot {input} | gzip > {output} 2> {log}
        """


BWA_INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]


rule bam_bwa_sra:
    """Map SRR data on the fly and convert to bam.

    NB: This will use the SampleName column to identify the sample
    (not the SRS id)."""
    output:
        cram="<project>/bam_bwa_sra/{samplealias}/{srrun}.sort.md.cram",
        crai="<project>/bam_bwa_sra/{samplealias}/{srrun}.sort.md.cram.crai",
    input:
        srr="<project>/sra/{srrun}/{srrun}.sra",
        reference=f"<project>/ref/{config['reference']}",
        index=expand(
            "<project>/{ref}.{sfx}",
            ref=f"ref/{config['reference']}",
            sfx=BWA_INDEX_SUFFIX,
        ),
    params:
        rg=get_read_group,
        options="-M",
        regions=config.get("regions", ""),
    conda:
        "../envs/bwamem.yaml"
    benchmark:
        "benchmarks/<project>/bam_bwa_sra/{samplealias}/{srrun}.sort.md.bam.benchmark.txt"
    log:
        "logs/<project>/bam_bwa_sra/{samplealias}/{srrun}.sort.md.bam.log",
    threads: 12
    priority: 50
    shell:
        """
        fasterq-dump --skip-technical -Z --split-spot {input.srr} |
            bwa mem {params.rg} -p -t {threads} {params.options} {input.reference} - 2>> {log}|
            samtools fixmate -m - /dev/stdout |
            samtools sort - | samtools markdup - /dev/stdout |
            samtools view -h -C -o {output.cram} > {log} 2>&1
        samtools index {output.cram} >> {log} 2>&1
        """
