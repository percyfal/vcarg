rule gatk_haplotypecaller:
    """Run GATK HaplotypeCaller"""
    output:
        vcf="<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz",
        tbi="<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz.tbi",
    input:
        bam=branch(
            lambda wildcards: wildcards.callmode == "raw",
            then="<project>/cram_2_bam/{samplealias}.bam",
            otherwise="<project>/gatk-bqsr/{samplealias}.bam",
        ),
        csi=branch(
            lambda wildcards: wildcards.callmode == "raw",
            then="<project>/cram_2_bam/{samplealias}.bam.csi",
            otherwise="<project>/gatk-bqsr/{samplealias}.bai",
        ),
        ref=f"<ref>/{config['reference']}",
        dict=re.sub("(.fasta|.fa)$", ".dict", f"<ref>/{config['reference']}"),
        fai=f"<ref>/{config['reference']}.fai",
    wildcard_constraints:
        mode="(.g|)",
    params:
        mode_options=lambda wildcards: "-ERC GVCF" if wildcards.mode == ".g" else "",
        options=" ".join(
            [
                "-A",
                "FisherStrand",
                "-A",
                "QualByDepth",
                "-A",
                "MappingQuality",
                "-G",
                "StandardAnnotation",
            ]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_haplotypecaller/<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_haplotypecaller/<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz.log",
    priority: 100
    threads: 4
    shell:
        """
        gatk HaplotypeCaller -OVI true \
        {params.options} {params.mode_options} \
        --input {input.bam} --output {output.vcf} \
        --reference {input.ref} \
        --native-pair-hmm-threads {threads} > {log} 2>&1
        """


rule gatk_raw_or_bqsr_variant_filtration:
    """Filter raw or bqsr variants"""
    output:
        vcf="<project>/gatk-hc-filter-{callmode}/{samplealias}{mode}.vcf.gz",
        tbi="<project>/gatk-hc-filter-{callmode}/{samplealias}{mode}.vcf.gz.tbi",
    input:
        vcf="<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz",
        tbi="<project>/gatk-hc-{callmode}/{samplealias}{mode}.vcf.gz.tbi",
    params:
        options=gatk_raw_or_bqsr_variant_filtration_options,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_raw_or_bqsr_variant_filtration/<project>/gatk-hc-filter-{callmode}/{samplealias}{mode}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_raw_or_bqsr_variant_filtration/<project>/gatk-hc-filter-{callmode}/{samplealias}{mode}.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk VariantFiltration -OVI true --variant {input.vcf} --output {output.vcf} {params.options} > {log} 2>&1
        """


rule gatk_base_recalibrator:
    """Recalibrate bases using raw variant calls as known sites"""
    output:
        table="<project>/gatk-bqsr/{samplealias}.table",
    input:
        bam="<project>/cram_2_bam/{samplealias}.bam",
        csi="<project>/cram_2_bam/{samplealias}.bam.csi",
        known_sites=expand(
            "<project>/gatk-hc-raw/{samplealias}.g.vcf.gz",
            samplealias=sampleinfo.SampleAlias.values,
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        known_sites=lambda wildcards, input: " ".join(
            [f"--known-sites {x}" for x in input.known_sites]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_base_recalibrator/<project>/gatk-bqsr/{samplealias}.table.benchmark.txt"
    log:
        "<logs>/gatk_base_recalibrator/<project>/gatk-bqsr/{samplealias}.table.log",
    threads: 1
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} {params.known_sites} -O {output.table} > {log} 2>&1
        """


rule gatk_apply_bqsr:
    """Apply BQSR on input bam"""
    output:
        recal=temp("<project>/gatk-bqsr/{samplealias}.bam"),
        bai=temp("<project>/gatk-bqsr/{samplealias}.bai"),
    input:
        table="<project>/gatk-bqsr/{samplealias}.table",
        bam="<project>/cram_2_bam/{samplealias}.bam",
        csi="<project>/cram_2_bam/{samplealias}.bam.csi",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_apply_bqsr/<project>/gatk-bqsr/{samplealias}.output.txt"
    log:
        "<logs>/gatk_apply_bqsr/<project>/gatk-bqsr/{samplealias}.log",
    threads: 1
    shell:
        """
        gatk ApplyBQSR -OBI -bqsr {input.table} -I {input.bam} -O {output.recal} > {log} 2>&1
        """


rule gatk_base_recalibrator_after:
    """Generate calibrated table for recalibrated BAM file"""
    output:
        table="<project>/gatk-bqsr/{samplealias}.after.table",
    input:
        bam="<project>/gatk-bqsr/{samplealias}.bam",
        bai="<project>/gatk-bqsr/{samplealias}.bai",
        known_sites=expand(
            "<project>/gatk-hc-filter-raw/{samplealias}.g.vcf.gz",
            samplealias=sampleinfo.SampleAlias.values,
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        known_sites=lambda wildcards, input: " ".join(
            [f"--known-sites {x}" for x in input.known_sites]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_base_recalibrator_after/<project>/gatk-bqsr/{samplealias}.after.table.benchmark.txt"
    log:
        "<logs>/gatk_base_recalibrator_after/<project>/gatk-bqsr/{samplealias}.after.table.log",
    threads: 1
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} {params.known_sites} -O {output.table} > {log} 2>&1
        """


rule gatk_analyze_covariates:
    """Analyze covariates from before and after bqsr"""
    output:
        csv="<project>/gatk-bqsr/{samplealias}.after.csv",
        pdf="<project>/gatk-bqsr/{samplealias}.after.pdf",
    input:
        before="<project>/gatk-bqsr/{samplealias}.table",
        after="<project>/gatk-bqsr/{samplealias}.after.table",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_analyze_covariates/<project>/gatk-bqsr/{samplealias}.after.csv.benchmark.txt"
    log:
        "<logs>/gatk_analyze_covariates/<project>/gatk-bqsr/{samplealias}.after.csv.log",
    threads: 1
    shell:
        """
        gatk AnalyzeCovariates --before-report-file {input.before} \
        --after-report-file {input.after} --plots-report-file {output.pdf} \
        --intermediate-csv-file {output.csv} > {log} 2>&1
        """


rule gatk_bqsr_bam_2_cram:
    """ConvertBQSR BAM to CRAM format"""
    output:
        cram="<project>/gatk-bqsr-cram/{samplealias}.cram",
        crai="<project>/gatk-bqsr-cram/{samplealias}.cram.crai",
    input:
        bam="<project>/gatk-bqsr/{samplealias}.bam",
        bai="<project>/gatk-bqsr/{samplealias}.bai",
        ref=f"<ref>/{config['reference']}",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "<benchmarks>/gatk_bqsr_bam_2_cram/<project>/gatk-bqsr/{samplealias}.cram.benchmark.txt"
    log:
        "<logs>/gatk_bqsr_bam_2_cram/<project>/gatk-bqsr/{samplealias}.cram.log",
    threads: 1
    shell:
        """
        samtools view --write-index -C -T {input.ref} -o {output.cram} {input.bam} > {log} 2>&1
        """


rule gatk_combine_gvcfs:
    """Run GATK CombineGVCFs.

    NB: the callmode=raw is currently not used but could be setup as
    input to BQSR as known sites.
    """
    output:
        vcf="<project>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz",
        tbi="<project>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.tbi",
    input:
        vcf=expand(
            "<project>/gatk-hc-{{callmode}}/{samplealias}.g.vcf.gz",
            samplealias=sampleinfo["SampleAlias"],
        ),
        tbi=expand(
            "<project>/gatk-hc-{{callmode}}/{samplealias}.g.vcf.gz.tbi",
            samplealias=sampleinfo["SampleAlias"],
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        vcf=lambda wildcards, input: " ".join([f"-V {x}" for x in input.vcf]),
        intervals=" ".join([f"-L {ivl}" for ivl in config.get("regions", [])]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_combine_gvcfs/<project>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_combine_gvcfs/<project>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk CombineGVCFs -OVI true {params.intervals} --output {output.vcf} --reference {input.ref} {params.vcf} > {log} 2>&1
        """


rule gatk_genotype_gvcfs:
    """GATK GenotypeGVCFs"""
    output:
        vcf="<project>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz",
        tbi="<project>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz.tbi",
    input:
        vcf="<project>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz",
        ref=f"<ref>/{config['reference']}",
    params:
        intervals=" ".join([f"-L {ivl}" for ivl in config.get("regions", [])]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_genotype_gvcfs/<project>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_genotype_gvcfs/<project>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.log",
    threads: 1
    shell:
        """
        gatk GenotypeGVCFs -OVI true {params.intervals} \
        -R {input.ref} -V {input.vcf} \
        -O {output.vcf} --all-sites \
        > {log} 2>&1
        """
