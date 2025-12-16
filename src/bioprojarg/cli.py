"""Command line interface for the application."""

import logging
import pathlib
import shutil
import subprocess
from importlib import resources

import click


from .env import Environment
from . import options


from . import __version__

logger = logging.getLogger(__name__)


def snakemake(*, targets=None, smk_options=None, snakefile=None):
    """Run snakemake workflows."""
    if isinstance(smk_options, list):
        smk_options = " ".join(smk_options)
    cmdlist = [
        "snakemake",
        f"{'-s ' + str(snakefile) if snakefile else ''}",
        f"{str(smk_options) or ''}",
        f"{str(targets) or ''}",
    ]
    cmd = " ".join(cmdlist)
    if shutil.which("snakemake") is None:
        logger.info("snakemake not installed; cannot run command:")
        logger.info("  %s", cmd)
        return

    try:
        logger.debug("running %s", cmd)
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError:
        logger.error("%s failed", cmd)
        raise


PKG_DIR = pathlib.Path(__file__).absolute().parent
CONTEXT_SETTINGS = {"auto_envvar_prefix": "RECRATE", "show_default": True}

pass_environment = click.make_pass_decorator(Environment, ensure=True)


@click.group(
    context_settings=CONTEXT_SETTINGS,
    help=__doc__,
    name="bioprojarg",
)
@click.version_option(version=__version__)
@options.debug_option()
@pass_environment
def cli(env):
    """CLI docstring for recrate"""
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s"
    )
    if env.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    logger.debug("Starting bioprojarg CLI")


@cli.command(
    context_settings={"ignore_unknown_options": True}, help="Run snakemake workflow"
)
@options.test_option()
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run(test, snakemake_args):
    """Run snakemake workflow"""
    click.echo("Running Snakemake workflow")
    smk_options = " ".join(list(snakemake_args) + test)
    snakefile = resources.files("bioprojarg") / "workflow" / "Snakefile"
    snakemake(smk_options=smk_options, snakefile=snakefile, targets="")
