rule download_datasources:
    """Download datasources listed in config file"""
    output:
        urltarget="<project>/{urltarget}",
    input:
        lambda wildcards: storage(
            config.get("datasources", {}).get(wildcards.urltarget)
        ),
    wildcard_constraints:
        urltarget=f'({"|".join([str(x) for x in set(datasources.keys())])})',
    conda:
        "../envs/storage.yaml"
    benchmark:
        "benchmarks/{urltarget}.benchmark.txt"
    log:
        "logs/{urltarget}.log",
    threads: 1
    shell:
        """
        mv {input} {output.urltarget}
        """
