#!/usr/bin/env python3
"""Tree sequence inference with tsinfer"""

from snakemake.script import snakemake

import logging
import random

import tsinfer  # pylint: disable=import-error
import zarr  # pylint: disable=import-error

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_tsinfer(
    zarrfile,
    outfile,
    chromosome,
    random_seed,
    recombination_rate,
    recombination_rate_map,
    mismatch_ratio,
    threads,
):  # pylint: disable=too-many-arguments, too-many-positional-arguments
    """Run tsinfer on variant data"""
    random.seed(random_seed)

    # Create chromosome mask
    vcz = zarr.open(zarrfile, mode="r")

    contig = vcz["contig_id"][:].tolist().index(chromosome)
    site_mask = vcz["variant_contig"][:] != contig

    variant_data = tsinfer.VariantData(
        zarrfile, ancestral_state="variant_AA", site_mask=site_mask
    )
    logger.info("There are %i variants in the zarr", variant_data.num_sites)

    anc_data = tsinfer.generate_ancestors(variant_data, progress_monitor=True)
    logger.info(
        "Generated %i ancestors from %i sites",
        anc_data.num_ancestors,
        variant_data.num_sites,
    )

    ancestors_arg = tsinfer.match_ancestors(
        variant_data, anc_data, progress_monitor=False
    )
    logger.info(ancestors_arg)

    if recombination_rate_map is not None:
        recombination_rate = recombination_rate_map
    arg = tsinfer.match_samples(
        variant_data,
        ancestors_arg,
        recombination_rate=recombination_rate,
        mismatch_ratio=mismatch_ratio,
        progress_monitor=True,
    )
    logger.info(arg)

    arg = tsinfer.infer(variant_data, num_threads=threads)
    arg.dump(outfile)


run_tsinfer(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.wildcards.chrom,
    snakemake.params.random_seed,
    snakemake.params.recombination_rate,
    snakemake.params.recombination_rate_map,
    snakemake.params.mismatch_ratio,
    snakemake.threads,
)
