import os
from os.path import dirname as dname
import pandas as pd
import numpy as np
from snakemake.utils import validate
from snakemake.logging import logger


# Validate config
validate(config, schema="../schemas/config.schema.yaml")

try:
    sampleinfo = pd.read_csv("sampleinfo.csv", sep=",")
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


class Interval:
    """Zero-based genomic interval."""

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    def __repr__(self):
        return f"{self.chrom}\t{self.start}\t{self.end}"

    @staticmethod
    def parse(input_str):
        fields = input_str.strip().split("\t")
        if len(fields) == 3:
            chrom, start, end = fields
        elif len(fields) == 1:
            chrom, positions = input_str.split(":")
            start, end = positions.split("-")
        else:
            chrom = fields[0]
            start = 0
            end = int(fields[1])
        return Interval(chrom, int(start), int(end))


def get_reference_basename():
    reference = config["reference"]
    reference_basename = re.sub(r"(\.fasta|\.fa|\.fna)(.gz|)$", "", reference)
    return reference_basename


def get_samplename_dict(wildcards):
    samplename = wildcards.samplename
    samplealias = sampleinfo.set_index("SampleName").loc[samplename].SampleAlias
    srrun = sampleinfo.set_index("SampleName").loc[samplename].Run
    info = {"srrun": srrun, "samplename": samplename, "samplealias": samplealias}
    return info


def get_read_group(wildcards):
    info = get_samplename_dict(wildcards)
    run = info["srrun"]
    sample = info["samplename"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_sample_runs(wildcards):
    samplename = wildcards.samplename
    srrun = sampleinfo.set_index("SampleName").loc[samplename].Run
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


def _make_intervals():
    """Split the reference genome into intervals based on config parameters."""
    logger.info(
        (
            "Creating intervals for variant calling and downstream "
            "processing from reference genome..."
        )
    )
    index_file = os.path.join("ref", config["reference"] + ".fai")
    if not os.path.exists(index_file):
        logger.warning(f"Reference index file {index_file} not found.")
        return []
    ivl = []
    min_interval_length = config.get("min_interval_length", 0)
    max_interval_length = config.get("max_interval_length", np.inf)
    n = {"total": 0, "passed": 0}
    if config.get("regions", None) is not None:
        logger.info("Using user-specified regions for interval creation.")
        source = config["regions"]
    else:
        source = open(index_file)
    for line in source:
        n["total"] = n["total"] + 1
        iv = Interval.parse(line)
        if len(iv) < min_interval_length:
            logger.warning(
                f"Skipping {iv.chrom} of length {len(iv)} < min_interval_length {min_interval_length}"
            )
            continue
        n["passed"] = n["passed"] + 1
        n_intervals = int(len(iv) / max_interval_length) + 1
        breaks = np.linspace(iv.start, iv.end, n_intervals + 1, dtype=int)
        for i in range(n_intervals):
            ivl.append(Interval(iv.chrom, breaks[i], breaks[i + 1]))
    logger.info(f"Kept {n['passed']} out of {n['total']} contigs")
    return ivl


def group_intervals():
    """Group intervals into batches based on interval batch size"""
    ivl = _make_intervals()
    if len(ivl) == 0:
        return {}
    n = {"total": len(ivl)}
    intervals = {}
    interval_batch_size = config.get("interval_batch_size", np.inf)
    cumsum = np.cumsum([len(reg) for reg in ivl])
    indices = np.array(cumsum / interval_batch_size, dtype=int)
    for j in range(max(indices) + 1):
        i = np.where(indices == j)
        group_id = f"{ivl[i[0][0]]}"
        intervals[group_id] = [ivl[k] for k in i[0]]
    logger.info(
        f"{n['total']} intervals grouped into {len(intervals)} chunks for parallel processing."
    )
    logger.info("Intervals:" + ",".join(intervals.keys()))
    return intervals
