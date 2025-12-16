rule download_sraruninfo:
    """Download sraruninfo for bioproject."""
    output:
        csv="<project>/resources/SraRunInfo.csv",
    params:
        bioproj=lambda wildcards, output: dname(dname(output.csv)),
        pyfilter=config.get("sraruninfo", {}).get("filter", ""),
    conda:
        "../envs/sratools.yaml"
    benchmark:
        "<benchmarks>/download_sraruninfo/<project>/resources/SraRunInfo.csv.benchmark.txt"
    log:
        "<logs>/download_sraruninfo/<project>/resources/SraRunInfo.csv.log",
    threads: 1
    shell:
        """
        esearch -db sra -query '{params.bioproj}' |
            efetch -format runinfo |
            qsv py filter "{params.pyfilter}" > {output.csv} 2> {log}
        """


rule make_samples_table:
    """Download sample information and create samples table"""
    output:
        csv="<project>/resources/samples.csv",
    input:
        storage.http(config.get("samplesheet", {}).get("url", [])),
    params:
        exe=config.get("samplesheet", {}).get("script", ""),
    conda:
        "../envs/textutils.yaml"
    benchmark:
        "<benchmarks>/make_samples_table/<project>/resources/samples.csv.benchmark.txt"
    log:
        "<logs>/make_samples_table/<project>/resources/samples.csv.log",
    threads: 1
    script:
        "{params.exe}"


rule join_sraruninfo_samples:
    """Join SraRunInfo.csv and samples.csv"""
    output:
        sampleinfo="<project>/sampleinfo.csv",
    input:
        runinfo="<project>/resources/SraRunInfo.csv",
        samples="<project>/resources/samples.csv",
    params:
        exe=config.get("sampleinfo_script", ""),
    conda:
        "../envs/textutils.yaml"
    benchmark:
        "benchmarks/sampleinfo.csv.benchmark.txt"
    log:
        "logs/sampleinfo.csv.log",
    threads: 1
    script:
        "{params.exe}"


localrules:
    make_samples_table,
    join_sraruninfo_samples,
