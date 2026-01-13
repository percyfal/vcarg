#!/usr/bin/env python3
"""Ancestral allele annotation with fastDFE"""

import random
import logging

import fastdfe as fd
from snakemake.script import snakemake
from snakemake.shell import shell


logger = logging.getLogger(__name__)


random.seed(snakemake.params.random_seed)

vcf = snakemake.input.vcf
output = snakemake.output.vcf
population_column = snakemake.params.population_column
sample_column = snakemake.params.sample_column
outgroup_populations = snakemake.params.outgroup_populations
n_ingroups = snakemake.params.n_ingroups
sampleinfo = snakemake.params.sampleinfo

outgroups = []
print(sampleinfo.head())
print(population_column)
print(sample_column)
print(outgroup_populations)

for pop in outgroup_populations:
    if pop not in sampleinfo[population_column].unique():
        raise ValueError(
            f"Outgroup population '{pop}' not found in samplesheet"
            f"column '{population_column}'"
        )
    outgroups.append(
        random.choice(
            sampleinfo[sampleinfo[population_column] == pop][sample_column].tolist()
        )
    )
ingroups = sampleinfo[~sampleinfo[population_column].isin(outgroup_populations)][
    sample_column
].tolist()
if len(ingroups) < n_ingroups:
    logger.warning(
        ("Not enough ingroup samples available (%i) to select %i ingroups",),
        len(ingroups),
        n_ingroups,
    )
    logger.info("Using all available ingroup samples as ingroups")

ann = fd.Annotator(
    vcf=vcf,
    annotations=[
        fd.MaximumLikelihoodAncestralAnnotation(
            outgroups=outgroups, n_ingroups=n_ingroups, ingroups=ingroups
        )
    ],
    output=output,
)

ann.annotate()

shell(f"bcftools index {snakemake.output.vcf}")
