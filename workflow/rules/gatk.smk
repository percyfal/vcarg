rule gatk_haplotypecaller_intervals:
    """Run GATK HaplotypeCaller on intervals"""
    output:
        vcf=temp("<work>/gatk-hc-{callmode}/{samplename}/{ivl}{mode}.vcf.gz"),
        tbi=temp("<work>/gatk-hc-{callmode}/{samplename}/{ivl}{mode}.vcf.gz.tbi"),
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
        seqdict=f"<ref>/{reference_basename}.dict",
        fai=f"<ref>/{config['reference']}.fai",
        intervals="<ref>/intervals/{ivl}.bed",
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
        "<benchmarks>/gatk_haplotypecaller/<work>/gatk-hc-{callmode}/{samplename}/{ivl}{mode}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_haplotypecaller/<work>/gatk-hc-{callmode}/{samplename}/{ivl}{mode}.vcf.gz.log",
    priority: 100
    threads: 4
    shell:
        """
        gatk HaplotypeCaller -OVI true \
        -L {input.intervals} \
        {params.options} {params.mode_options} \
        --input {input.bam} --output {output.vcf} \
        --reference {input.ref} \
        --native-pair-hmm-threads {threads} > {log} 2>&1
        """


rule gatk_raw_or_bqsr_variant_filtration:
    """Filter raw or bqsr variants"""
    output:
        vcf=temp("<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz"),
        tbi=temp("<work>/gatk-hc-filter-{callmode}/{samplename}{mode}.vcf.gz.tbi"),
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
        known_sites="<results>/gatk-select-filtered-variants-raw/all.vcf.gz",
        known_sites_tbi="<results>/gatk-select-filtered-variants-raw/all.vcf.gz.tbi",
        ref=f"<ref>/{config['reference']}",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_base_recalibrator/<work>/gatk-bqsr/{samplename}.table.benchmark.txt"
    log:
        "<logs>/gatk_base_recalibrator/<work>/gatk-bqsr/{samplename}.table.log",
    threads: 2
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output.table} > {log} 2>&1
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
    threads: 2
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
        known_sites="<results>/gatk-genotype-gvcf-filter-raw/all.vcf.gz",
        known_sites_tbi="<results>/gatk-genotype-gvcf-filter-raw/all.vcf.gz.tbi",
        # known_sites=expand(
        #     "<work>/gatk-hc-filter-raw/{samplename}.g.vcf.gz",
        #     samplename=sampleinfo.SampleName.values,
        # ),
        # known_sites_tbi=expand(
        #     "<work>/gatk-hc-filter-raw/{samplename}.g.vcf.gz.tbi",
        #     samplename=sampleinfo.SampleName.values,
        # ),
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
    threads: 2
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
    threads: 4
    shell:
        """
        gatk AnalyzeCovariates --before-report-file {input.before} \
        --after-report-file {input.after} --plots-report-file {output.pdf} \
        --intermediate-csv-file {output.csv} > {log} 2>&1
        """


rule gatk_bqsr_bam_2_cram:
    """Convert BQSR BAM to CRAM format"""
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


rule gatk_sample_name_map:
    """Map sample names to file paths for GenomicsDB import."""
    output:
        txt="<work>/gatk-sample-name-map/{ivl}{mode}-{callmode}.samples.txt",
    input:
        bed="<ref>/intervals/{ivl}.bed",
        vcf=expand(
            "<work>/gatk-hc-{{callmode}}/{samplename}/{{ivl}}{{mode}}.vcf.gz",
            samplename=sampleinfo.SampleName.values,
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_sample_name_map/<work>/gatk-sample-name-map/{ivl}{mode}-{callmode}.samples.txt.benchmark.txt"
    log:
        "<logs>/gatk_sample_name_map/<work>/gatk-sample-name-map/{ivl}{mode}-{callmode}.samples.txt.log",
    threads: 1
    shell:
        """
        for vcf in {input.vcf}; do
            sample=$(basename $(dirname $vcf))
            echo -e "${{sample}}\t${{vcf}}"
        done > {output.txt}
        """


rule gatk_genomics_db_import_intervals:
    """Import intervals into GenomicsDB for GATK joint genotyping."""
    output:
        directory("<work>/gatk-genomics-db-import-intervals-{callmode}/{ivl}{mode}"),
    input:
        intervals="<ref>/intervals/{ivl}.bed",
        sample_name_map="<work>/gatk-sample-name-map/{ivl}{mode}-{callmode}.samples.txt",
        vcf=expand(
            "<work>/gatk-hc-{{callmode}}/{samplename}/{{ivl}}{{mode}}.vcf.gz",
            samplename=sampleinfo.SampleName.values,
        ),
        tbi=expand(
            "<work>/gatk-hc-{{callmode}}/{samplename}/{{ivl}}{{mode}}.vcf.gz.tbi",
            samplename=sampleinfo.SampleName.values,
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_genomics_db_import_intervals/<work>/gatk-genomics-db-import-intervals-{callmode}/{ivl}{mode}.benchmark.txt"
    log:
        "<logs>/gatk_genomics_db_import_intervals/<work>/gatk-genomics-db-import-intervals-{callmode}/{ivl}{mode}.log",
    threads: 14
    priority: 200
    shell:
        """
        gatk GenomicsDBImport -OVI true --genomicsdb-workspace-path {output} \
        -L {input.intervals} \
        --sample-name-map {input.sample_name_map} \
        --batch-size 50 \
        --reader-threads {threads} > {log} 2>&1
        """


rule gatk_combine_gvcfs:
    """Run GATK CombineGVCFs.

    NB: the callmode=raw is currently not used but could be setup as
    input to BQSR as known sites.
    """
    output:
        vcf=temp("<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz"),
        tbi=temp("<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz.tbi"),
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


rule gatk_genotype_gvcfs_intervals:
    """GATK GenotypeGVCFs by intervals using GenomicsDB as input."""
    output:
        vcf=temp(
            "<work>/gatk-genotype-gvcf-{callmode}/{callset}{allsites}.{ivl}.vcf.gz"
        ),
        tbi=temp(
            "<work>/gatk-genotype-gvcf-{callmode}/{callset}{allsites}.{ivl}.vcf.gz.tbi"
        ),
    input:
        # FIXME: alternative input would be from combine-GVCFs:
        # vcf="<work>/gatk-combine-gvcf-{callmode}/{callset}.g.vcf.gz"
        db="<work>/gatk-genomics-db-import-intervals-{callmode}/{ivl}.g",
        ref=f"<ref>/{config['reference']}",
        intervals="<ref>/intervals/{ivl}.bed",
    params:
        allsites=lambda wildcards: (
            "--include-non-variant-sites" if wildcards.allsites == ".allsites" else ""
        ),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_genotype_gvcfs/<work>/gatk-genotype-gvcf-{callmode}/{callset}{allsites}{ivl}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_genotype_gvcfs/<work>/gatk-genotype-gvcf-{callmode}/{callset}{allsites}{ivl}.log",
    threads: 1
    shell:
        """
        gatk GenotypeGVCFs -OVI true -L {input.intervals} \
        -R {input.ref} -V gendb://{input.db} \
        -O {output.vcf} {params.allsites} \
        > {log} 2>&1
        """


rule gatk_gather_vcfs:
    """Gather GATK GenotypeGVCFs interval VCFs into single VCF per callset."""
    output:
        vcf="<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz",
    input:
        vcf=expand(
            "<work>/gatk-genotype-gvcf-{{callmode}}/{{callset}}{{allsites}}.{ivl}.vcf.gz",
            ivl=intervals.keys(),
        ),
        tbi=expand(
            "<work>/gatk-genotype-gvcf-{{callmode}}/{{callset}}{{allsites}}.{ivl}.vcf.gz.tbi",
            ivl=intervals.keys(),
        ),
    params:
        vcf=lambda wildcards, input: " ".join([f"--INPUT {x}" for x in input.vcf]),
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_gather_vcfs/<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_gather_vcfs/<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk GatherVcfs {params.vcf} --OUTPUT {output.vcf} > {log} 2>&1
        """


rule bcftools_create_index:
    """Create index with bcftools"""
    output:
        tbi="<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.tbi",
    input:
        vcf="<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "<benchmarks>/bcftools_create_index/<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.tbi.benchmark.txt"
    log:
        "<logs>/bcftools_create_index/<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.tbi.log",
    threads: 1
    shell:
        """
        bcftools index -t {input.vcf}
        """


rule gatk_raw_or_bqsr_genotype_variant_filtration:
    """Annotate raw or bqsr variants from GATK genotyped GVCFs with filter information."""
    output:
        vcf="<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz",
        tbi="<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz.tbi",
    input:
        vcf="<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz",
        tbi="<results>/gatk-gather-vcfs-{callmode}/{callset}{allsites}.vcf.gz.tbi",
    params:
        options=gatk_raw_or_bqsr_variant_filtration_options,
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_raw_or_bqsr_genotype_variant_filtration/<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_raw_or_bqsr_genotype_variant_filtration/<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk VariantFiltration -OVI true --variant {input.vcf} --output {output.vcf} {params.options} > {log} 2>&1
        """


rule gatk_select_filtered_variants:
    """Select filtered variants from GATK genotyped and filtered VCFs"""
    output:
        vcf="<results>/gatk-select-filtered-variants-{callmode}/{callset}{allsites}.vcf.gz",
        tbi="<results>/gatk-select-filtered-variants-{callmode}/{callset}{allsites}.vcf.gz.tbi",
    input:
        vcf="<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz",
        tbi="<results>/gatk-genotype-gvcf-filter-{callmode}/{callset}{allsites}.vcf.gz.tbi",
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "<benchmarks>/gatk_select_filtered_variants/<results>/gatk-select-filtered-variants-{callmode}/{callset}{allsites}.vcf.gz.benchmark.txt"
    log:
        "<logs>/gatk_select_filtered_variants/<results>/gatk-select-filtered-variants-{callmode}/{callset}{allsites}.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk SelectVariants -OVI true --variant {input.vcf} --output {output.vcf} --exclude-filtered > {log} 2>&1
        """
