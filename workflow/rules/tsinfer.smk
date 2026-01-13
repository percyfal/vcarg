rule make_sampleinfo_population:
    """Make sampleinfo and population files for tsinfer"""
    output:
        sampleinfo="sampleinfo.tsinfer.csv",
        populations="populations.tsinfer.csv",
    input:
        csv="sampleinfo.csv",
    params:
        sampleinfo_columns=config.get("tsinfer", {}).get(
            "sampleinfo_columns", "Sample,Population"
        ),
        populations_columns=config.get("tsinfer", {}).get(
            "populations_columns", "Population"
        ),
        popindex=config.get("tsinfer", {})
        .get("populations_columns", "Population")
        .split(",")[0],
    conda:
        "../envs/textutils.yaml"
    benchmark:
        "<benchmarks>/make_sampleinfo_population/sampleinfo.tsinfer.csv.benchmark.txt"
    log:
        "<logs>/make_sampleinfo_population/sampleinfo.tsinfer.csv.log",
    threads: 1
    shell:
        """
        csvtk cut -f {params.sampleinfo_columns} {input.csv} > {output.sampleinfo} 2> {log}
        csvtk cut -f {params.populations_columns} {input.csv} | csvtk sort -k {params.popindex} | csvtk uniq -f {params.popindex} > {output.populations} 2>> {log}
        """


rule vcf_to_vcz:
    """Convert VCF to VCZ format"""
    output:
        vcz=directory("results/vcz/{callset}.{chrom}.vcz"),
    input:
        vcf="<results>/ancestral_allele/{callset}.vcf.gz",
        csi="<results>/ancestral_allele/{callset}.vcf.gz.csi",
    conda:
        "../envs/bio2zarr.yaml"
    benchmark:
        "<benchmarks>/vcf_to_csv/results/vcz/{callset}.{chrom}.vcz.benchmark.txt"
    log:
        "<logs>/vcf_to_csv/results/vcz/{callset}.{chrom}.vcz.log",
    threads: 1
    shell:
        """
        python -m bio2zarr vcf2zarr convert -p {threads} --force {input.vcf} {output.vcz} > {log} 2>&1
        """


rule tsinfer_inference:
    """Run tsinfer to infer tree sequence from phased VCF"""
    output:
        trees="<results>/tsinfer/{callset}.{chrom}.trees",
    input:
        vcz="<results>/vcz/{callset}.{chrom}.vcz",
        sampleinfo="sampleinfo.tsinfer.csv",
        populations="populations.tsinfer.csv",
    params:
        random_seed=42,
        mutation_rate=config.get("mutation_rate", 1.25e-8),
        recombination_rate=config.get("recombination_rate", 1e-8),
        recombination_rate_map=config.get("recombination_rate_map", None),
        mismatch_ratio=config.get("tsinfer", {}).get("mismatch_ratio", 0.001),
    conda:
        "../envs/tsinfer.yaml"
    benchmark:
        "<benchmarks>/tsinfer_inference/<results>/tsinfer/{callset}.{chrom}.trees.benchmark.txt"
    log:
        "<logs>/tsinfer_inference/<results>/tsinfer/{callset}.{chrom}.trees.log",
    threads: 15
    script:
        "../scripts/tsinfer_infer.py"


rule tsdate:
    """Run tsdate to date tree sequence"""
    output:
        trees="<results>/tsinfer/{callset}.{chrom}.dated.trees",
    input:
        trees="<results>/tsinfer/{callset}.{chrom}.preprocessed.trees",
    params:
        mutation_rate=config.get("mutation_rate", 1.25e-8),
        recombination_rate=config.get("recombination_rate", 1e-8),
    conda:
        "../envs/tsinfer.yaml"
    benchmark:
        "<benchmarks>/tsdate/<results>/tsinfer/{callset}.{chrom}.dated.trees.benchmark.txt"
    log:
        "<logs>/tsdate/<results>/tsinfer/{callset}.{chrom}.dated.trees.log",
    threads: 1
    shell:
        """
        tsdate date -m {params.mutation_rate} -p {input.trees} {output.trees} > {log} 2>&1
        """


rule tspreprocess:
    """Preprocess tree sequence for analysis with tskit"""
    output:
        trees="<results>/tsinfer/{callset}.{chrom}.preprocessed.trees",
    input:
        trees="<results>/tsinfer/{callset}.{chrom}.trees",
    conda:
        "../envs/tsinfer.yaml"
    benchmark:
        "<benchmarks>/tspreprocess/<results>/tsinfer/{callset}.{chrom}.preprocessed.trees.benchmark.txt"
    log:
        "<logs>/tspreprocess/<results>/tsinfer/{callset}.{chrom}.preprocessed.trees.log",
    threads: 1
    shell:
        """
        tsdate preprocess {input.trees} {output.trees} > {log} 2>&1
        """
