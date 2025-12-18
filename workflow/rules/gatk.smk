rule gatk_haplotypecaller:
    """Run GATK HaplotypeCaller"""
    output:
        vcf="<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz",
        tbi="<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz.tbi",
    input:
        bam=branch(
            lambda wildcards: wildcards.callmode == "raw",
            then="<work>/cram_2_bam/{samplename}.bam",
            otherwise="<work>/gatk-bqsr/{samplename}.bam",
        ),
        csi=branch(
            lambda wildcards: wildcards.callmode == "raw",
            then="<work>/cram_2_bam/{samplename}.bam.csi",
            otherwise="<work>/gatk-bqsr/{samplename}.bai",
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
        "<benchmarks>/gatk_haplotypecaller/<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_haplotypecaller/<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz.log",
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
        vcf="<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz",
        tbi="<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz.tbi",
    input:
        vcf="<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz",
        tbi="<work>/gatk-hc-{callmode}/{samplename}{mode}.vcf.gz.tbi",
    params:
        options=gatk_raw_or_bqsr_variant_filtration_options,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_raw_or_bqsr_variant_filtration/<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_raw_or_bqsr_variant_filtration/<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk VariantFiltration -OVI true --variant {input.vcf} --output {output.vcf} {params.options} > {log} 2>&1
        """


rule gatk_base_recalibrator:
    """Recalibrate bases using raw variant calls as known sites"""
    output:
        table="<work>/gatk-bqsr/{samplename}.table",
    input:
        bam="<work>/cram_2_bam/{samplename}.bam",
        csi="<work>/cram_2_bam/{samplename}.bam.csi",
        known_sites=expand(
            "<work>/gatk-hc-raw/{samplename}.g.vcf.gz",
            samplename=sampleinfo.SampleName.values,
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        known_sites=lambda wildcards, input: " ".join(
            [f"--known-sites {x}" for x in input.known_sites]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_base_recalibrator/<work>/gatk-bqsr/{samplename}.table.benchmark.txt"
    log:
        "<logs>/gatk_base_recalibrator/<work>/gatk-bqsr/{samplename}.table.log",
    threads: 1
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} {params.known_sites} -O {output.table} > {log} 2>&1
        """


rule gatk_apply_bqsr:
    """Apply BQSR on input bam"""
    output:
        recal=temp("<work>/gatk-bqsr/{samplename}.bam"),
        bai=temp("<work>/gatk-bqsr/{samplename}.bai"),
    input:
        table="<work>/gatk-bqsr/{samplename}.table",
        bam="<work>/cram_2_bam/{samplename}.bam",
        csi="<work>/cram_2_bam/{samplename}.bam.csi",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_apply_bqsr/<work>/gatk-bqsr/{samplename}.output.txt"
    log:
        "<logs>/gatk_apply_bqsr/<work>/gatk-bqsr/{samplename}.log",
    threads: 1
    shell:
        """
        gatk ApplyBQSR -OBI -bqsr {input.table} -I {input.bam} -O {output.recal} > {log} 2>&1
        """


rule gatk_base_recalibrator_after:
    """Generate calibrated table for recalibrated BAM file"""
    output:
        table="<work>/gatk-bqsr/{samplename}.after.table",
    input:
        bam="<work>/gatk-bqsr/{samplename}.bam",
        bai="<work>/gatk-bqsr/{samplename}.bai",
        known_sites=expand(
            "<work>/gatk-hc-filter-raw/{samplename}.g.vcf.gz",
            samplename=sampleinfo.SampleName.values,
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        known_sites=lambda wildcards, input: " ".join(
            [f"--known-sites {x}" for x in input.known_sites]
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_base_recalibrator_after/<work>/gatk-bqsr/{samplename}.after.table.benchmark.txt"
    log:
        "<logs>/gatk_base_recalibrator_after/<work>/gatk-bqsr/{samplename}.after.table.log",
    threads: 1
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} {params.known_sites} -O {output.table} > {log} 2>&1
        """


rule gatk_analyze_covariates:
    """Analyze covariates from before and after bqsr"""
    output:
        csv="<work>/gatk-bqsr/{samplename}.after.csv",
        pdf="<work>/gatk-bqsr/{samplename}.after.pdf",
    input:
        before="<work>/gatk-bqsr/{samplename}.table",
        after="<work>/gatk-bqsr/{samplename}.after.table",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_analyze_covariates/<work>/gatk-bqsr/{samplename}.after.csv.benchmark.txt"
    log:
        "<logs>/gatk_analyze_covariates/<work>/gatk-bqsr/{samplename}.after.csv.log",
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
        cram="<work>/gatk-bqsr-cram/{samplename}.cram",
        crai="<work>/gatk-bqsr-cram/{samplename}.cram.crai",
    input:
        bam="<work>/gatk-bqsr/{samplename}.bam",
        bai="<work>/gatk-bqsr/{samplename}.bai",
        ref=f"<ref>/{config['reference']}",
    conda:
        "../envs/bwa.yaml"
    benchmark:
        "<benchmarks>/gatk_bqsr_bam_2_cram/<work>/gatk-bqsr/{samplename}.cram.benchmark.txt"
    log:
        "<logs>/gatk_bqsr_bam_2_cram/<work>/gatk-bqsr/{samplename}.cram.log",
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
        vcf="<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz",
        tbi="<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.tbi",
    input:
        vcf=expand(
            "<work>/gatk-hc-{{callmode}}/{samplename}.g.vcf.gz",
            samplename=sampleinfo["SampleName"],
        ),
        tbi=expand(
            "<work>/gatk-hc-{{callmode}}/{samplename}.g.vcf.gz.tbi",
            samplename=sampleinfo["SampleName"],
        ),
        ref=f"<ref>/{config['reference']}",
    params:
        vcf=lambda wildcards, input: " ".join([f"-V {x}" for x in input.vcf]),
        intervals=" ".join([f"-L {ivl}" for ivl in config.get("regions", [])]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_combine_gvcfs/<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_combine_gvcfs/<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk CombineGVCFs -OVI true {params.intervals} --output {output.vcf} --reference {input.ref} {params.vcf} > {log} 2>&1
        """


rule gatk_genotype_gvcfs:
    """GATK GenotypeGVCFs"""
    output:
        vcf="<results>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz",
        tbi="<results>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz.tbi",
    input:
        vcf="<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz",
        ref=f"<ref>/{config['reference']}",
    params:
        intervals=" ".join([f"-L {ivl}" for ivl in config.get("regions", [])]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_genotype_gvcfs/<results>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_genotype_gvcfs/<results>/gatk-genotype-gvcf-{callmode}/{callset}.allsites.log",
    threads: 1
    shell:
        """
        gatk GenotypeGVCFs -OVI true {params.intervals} \
        -R {input.ref} -V {input.vcf} \
        -O {output.vcf} --all-sites \
        > {log} 2>&1
        """
