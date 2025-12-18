rule beagle_phase:
    """Phase and impute genotypes with beagle"""
    output:
        vcf="<results>/phasing/{callset}.vcf.gz",
        csi="<results>/phasing/{callset}.vcf.gz.csi",
    input:
        "<results>/biallelic-bqsr/{callset}.vcf.gz",
    resources:
        mem_mb=60000,
    params:
        prefix=lambda wildcards, output: output[0].replace(".vcf.gz", ""),
    conda:
        "../envs/phasing.yaml"
    benchmark:
        "<benchmarks>/beagle_phase/<results>/phasing/{callset}.vcf.gz.benchmark.txt"
    log:
        "<logs>/beagle_phase/<results>/phasing/{callset}.vcf.gz.log",
    threads: 15
    shell:
        """
        beagle -Xms{resources.mem_mb}m -Xmx{resources.mem_mb}m gt={input} out={params.prefix} nthreads={threads} > {log} 2>&1
        bcftools index {output.vcf}
        """


rule fastdfe_infer_aa:
    """Infer ancestral allele with FastDFE"""
    output:
        vcf="<results>/ancestral_allele/{callset}.vcf.gz",
        csi="<results>/ancestral_allele/{callset}.vcf.gz.csi",
    input:
        vcf="<results>/phasing/{callset}.vcf.gz",
        csi="<results>/phasing/{callset}.vcf.gz.csi",
    conda:
        "../envs/fastdfe.yaml"
    params:
        outgroups=" ".join(
            [f"--outgroup-populations {pop}" for pop in config.get("outgroups", [])]
        ),
        random_seed=config.get("fastdfe_random_seed", 42),
        population_column=config.get("population_column", "Population"),
        sample_column=config.get("sample_column", "SampleName"),
        outgroup_populations=config.get("outgroup_populations", []),
        n_ingroups=config.get("n_ingroups", 10),
        sampleinfo=sampleinfo,
    benchmark:
        "<benchmarks>/fastdfe_infer_aa/<results>/ancestral_allele/{callset}.vcf.gz.benchmark.txt"
    log:
        "<logs>/fastdfe_infer_aa/<results>/ancestral_allele/{callset}.vcf.gz.log",
    threads: 1
    script:
        "../scripts/fastdfe_aa.py"
