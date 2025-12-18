rule download_sraruninfo:
    """Download sraruninfo for bioproject."""
    output:
        csv="<resources>/SraRunInfo.<bioproject>.csv",
    params:
        bioproj=lambda wildcards, output: Path(output.csv).stem.split(".")[1],
        pyfilter=config.get("sraruninfo", {}).get("filter", ""),
    conda:
        "../envs/sratools.yaml"
    benchmark:
        "<benchmarks>/download_sraruninfo/<resources>/SraRunInfo.<bioproject>.csv.benchmark.txt"
    log:
        "<logs>/download_sraruninfo/<resources>/SraRunInfo.<bioproject>.csv.log",
    priority: 500
    threads: 1
    shell:
        """
        esearch -db sra -query '{params.bioproj}' |
            efetch -format runinfo > {output.csv}.bak 2> {log}
        cat {output.csv}.bak | qsv py filter "{params.pyfilter}" > {output.csv} 2>> {log}
        rm {output.csv}.bak
        """


rule make_samples_table:
    """Download sample information and create samples table"""
    output:
        csv="<resources>/samples.csv",
    input:
        storage.http(config.get("samplesheet", {}).get("url", [])),
    params:
        exe=config.get("samplesheet", {}).get("script", ""),
    conda:
        "../envs/textutils.yaml"
    benchmark:
        "<benchmarks>/make_samples_table/<resources>/samples.csv.benchmark.txt"
    log:
        "<logs>/make_samples_table/<resources>/samples.csv.log",
    threads: 1
    script:
        "{params.exe}"


rule join_sraruninfo_samples:
    """Join SraRunInfo.csv and samples.csv"""
    output:
        sampleinfo="sampleinfo.csv",
    input:
        runinfo="<resources>/SraRunInfo.<bioproject>.csv",
        samples="<resources>/samples.csv",
    params:
        exe=config.get("sampleinfo_script", ""),
    conda:
        "../envs/textutils.yaml"
    benchmark:
        "benchmarks/sampleinfo.csv.benchmark.txt"
    log:
        "logs/sampleinfo.csv.log",
    priority: 1000
    threads: 1
    script:
        "{params.exe}"


localrules:
    make_samples_table,
    join_sraruninfo_samples,
