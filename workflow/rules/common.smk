import os
from os.path import dirname as dname
import pandas as pd
from snakemake.utils import validate
from snakemake.logging import logger


# Validate config
validate(config, schema="../schemas/config.schema.yaml")

try:
    sampleinfo = pd.read_csv(bioproject / "sampleinfo.csv", sep=",")
except FileNotFoundError:
    logger.warning("sampleinfo.csv not found. Creating empty sampleinfo DataFrame.")
    logger.warning("Rerun the workflow once sampleinfo.csv has been generated")
    sampleinfo = pd.DataFrame(
        columns=[
            "SampleAlias",
            "SampleName",
            "Run",
        ]
    )


datasources = dict()
if "datasources" in config.keys():
    datasources = dict(
        zip(config["datasources"].keys(), config["datasources"].values())
    )


def get_samplealias_dict(wildcards):
    samplealias = wildcards.samplealias
    samplename = sampleinfo.set_index("SampleAlias").loc[samplealias].SampleName
    srrun = sampleinfo.set_index("SampleAlias").loc[samplealias].Run
    info = {"srrun": srrun, "samplename": samplename, "sample": samplealias}
    return info


def get_read_group(wildcards):
    info = get_samplealias_dict(wildcards)
    run = info["srrun"]
    sample = info["samplename"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_sample_runs(wildcards):
    samplealias = wildcards.samplealias
    srrun = sampleinfo.set_index("SampleAlias").loc[samplealias].Run
    return srrun


def gatk_raw_or_bqsr_variant_filtration_options(wildcards):
    if wildcards.callmode == "raw":
        options = [
            "--filter-name",
            "FisherStrand",
            "--filter",
            "'FS > 50.0'",
            "--filter-name",
            "QualByDepth",
            "--filter",
            "'QD < 4.0'",
            "--filter-name",
            "MappingQuality",
            "--filter",
            "'MQ < 50.0'",
        ]
    else:
        options = [
            "--filter-name",
            "FisherStrand",
            "--filter",
            "'FS > 60.0'",
            "--filter-name",
            "QualByDepth",
            "--filter",
            "'QD < 2.0'",
            "--filter-name",
            "MappingQuality",
            "--filter",
            "'MQ < 40.0'",
        ]

    return " ".join(options)
