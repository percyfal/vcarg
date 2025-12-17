BWA_INDEX_SUFFIX = ["amb", "ann", "bwt", "pac", "sa"]


rule bwa_index:
    """Make bwa index"""
    output:
        index=expand("<project>/bwa_index/{{prefix}}.fasta.{sfx}", sfx=BWA_INDEX_SUFFIX),
    input:
        fasta="<ref>/{prefix}.fasta",
    conda:
        "../envs/bwa.yaml"
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.index[0])[0],
    benchmark:
        "benchmarks/<project>/bwa_index/{prefix}.fasta.benchmark.txt"
    log:
        "logs/<project>/bwa_index/{prefix}.fasta.log",
    threads: 1
    shell:
        """
        bwa index -p {params.prefix} {input.fasta} > {log} 2>&1
        """


rule bam_bwa_srr:
    """Map SRR data on the fly and convert to bam, without marking duplicates.

    NB: This will use the SampleName column to identify the sample
    (not the SRS id)."""
    output:
        cram=temp("<project>/bam_bwa_srr/{samplealias}/{srrun}.sort.cram"),
        crai=temp("<project>/bam_bwa_srr/{samplealias}/{srrun}.sort.cram.crai"),
    input:
        srr="<project>/sra/{srrun}/{srrun}.sra",
        reference=f"<ref>/{config['reference']}",
        index=expand(
            "<project>/{ref}.{sfx}",
            ref=f"bwa_index/{config['reference']}",
            sfx=BWA_INDEX_SUFFIX,
        ),
    params:
        rg=get_read_group,
        options="-M",
        regions=config.get("regions", ""),
        index=lambda wildcards, input: os.path.splitext(input.index[0])[0],
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "benchmarks/<project>/bam_bwa_srr/{samplealias}/{srrun}.sort.cram.benchmark.txt"
    log:
        "logs/<project>/bam_bwa_srr/{samplealias}/{srrun}.sort.cram.log",
    threads: 12
    priority: 50
    shell:
        """
        fasterq-dump -e {threads} --skip-technical -Z --split-spot {input.srr} |
            bwa mem {params.rg} -p -t {threads} {params.options} {params.index} - 2> {log}|
            samtools fixmate -m - /dev/stdout |
            samtools sort - -T {resources.tmpdir} |
            samtools view -T {input.reference} -h -C -o {output.cram} >> {log} 2>&1
        samtools index {output.cram} >> {log} 2>&1
        if [[ "{params.regions}" != "" ]]; then
            samtools view -T {input.reference} -h -C -o {output.cram}.tmp {output.cram} {params.regions} >> {log} 2>&1
            mv -f {output.cram}.tmp {output.cram}
            samtools index {output.cram} >> {log} 2>&1
        fi
        """


rule merge_cram:
    """Merge run mappings to sample level cram"""
    output:
        cram="<project>/merge_cram/{samplealias}.cram",
        crai="<project>/merge_cram/{samplealias}.cram.crai",
    input:
        cram=lambda wildcards: expand(
            "<project>/bam_bwa_srr/{{samplealias}}/{srrun}.sort.cram",
            srrun=get_sample_runs(wildcards),
        ),
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "<benchmarks>/merge_cram/<project>/merge_cram/{samplealias}.cram.benchmark.txt"
    log:
        "<logs>/merge_cram/<project>/merge_cram/{samplealias}.cram.log",
    threads: 1
    shell:
        """
        samtools merge {input.cram} -o {output.cram} -O CRAM --write-index > {log} 2>&1
        """


rule cram_2_bam:
    """Convert CRAM to BAM format.

    Motivation: latest samtools outputs CRAM version 3.1 which is
    currently not supported by GATK.
    """
    output:
        bam=temp("<project>/cram_2_bam/{samplealias}.bam"),
        csi=temp("<project>/cram_2_bam/{samplealias}.bam.csi"),
    input:
        cram="<project>/merge_cram/{samplealias}.cram",
        crai="<project>/merge_cram/{samplealias}.cram.crai",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "<benchmarks>/<project>/cram_2_bam/{samplealias}.benchmark.txt"
    log:
        "<logs>/<project>/cram_2_bam/{samplealias}.log",
    threads: 1
    shell:
        """
        samtools view {input.cram} --write-index -O BAM -o {output.bam} > {log} 2>&1
    """
