"""Decorators for CLI"""

import logging
import typing as t

import click
from click.core import Command, Context, Parameter
from click.decorators import option

from vcarg.env import Environment

logger = logging.getLogger(__name__)

F = t.TypeVar("F", bound=t.Callable[..., t.Any])
FC = t.TypeVar("FC", t.Callable[..., t.Any], Command)


def debug_option(*param_decls: str, **kwargs: t.Any) -> t.Callable[[FC], FC]:
    """Add a ``--debug`` option which turns on debugging.

    :param param_decls: One or more option names. Defaults to the single
        value ``"--debug"``.
    :param kwargs: Extra arguments are passed to :func:`option`.
    """

    def callback(
        ctx: Context,
        param: Parameter,  # pylint: disable=unused-argument
        value: bool,  # pylint: disable=unused-argument
    ) -> None:
        if not value or ctx.resilient_parsing:
            return
        ctx.ensure_object(Environment)
        ctx.obj.debug = value
        if ctx.resilient_parsing:
            return

    if not param_decls:
        param_decls = ("--debug",)

    kwargs.setdefault("is_flag", True)
    kwargs.setdefault("expose_value", False)
    kwargs.setdefault("is_eager", True)
    kwargs.setdefault("help", ("Print debugging information."))
    kwargs["callback"] = callback
    return option(*param_decls, **kwargs)


def cores_option(default=None) -> t.Callable[[FC], FC]:
    """Add cores option."""

    def cores_callback(
        ctx: click.core.Context,  # pylint: disable=unused-argument
        param: click.core.Option,  # pylint: disable=unused-argument
        value: int,
    ) -> int:
        """Cores callback."""
        if value is None:
            return []
        if value < 1:
            logging.error("Cores must be greater than 0")
            raise ValueError("Cores must be greater than 0")
        return ["--cores", str(value)]

    return click.option(
        "-c",
        "--cores",
        help="number of cores",
        default=default,
        callback=cores_callback,
        type=int,
    )


def jobs_option(default=None) -> t.Callable[[FC], FC]:
    """Add jobs option."""

    def jobs_callback(
        ctx: click.core.Context,  # pylint: disable=unused-argument
        param: click.core.Option,  # pylint: disable=unused-argument
        value: int,
    ) -> int:
        """Jobs callback."""
        if value is None:
            return []
        if value < 1:
            logging.error("Jobs must be greater than 0")
            raise ValueError("Jobs must be greater than 0")
        return ["--jobs", str(value)]

    return click.option(
        "-j",
        "--jobs",
        help="number of jobs",
        default=default,
        callback=jobs_callback,
        type=int,
    )


def threads_option(default=None) -> t.Callable[[FC], FC]:
    """Add threads option."""

    def threads_callback(
        ctx: click.core.Context,  # pylint: disable=unused-argument
        param: click.core.Option,  # pylint: disable=unused-argument
        value: int,
    ) -> int:
        """Threads callback."""
        if value is None:
            return []
        if value < 1:
            logging.error("Threads must be greater than 0")
            raise ValueError("Threads must be greater than 0")
        return ["--threads", str(value)]

    return click.option(
        "-t",
        "--threads",
        help="number of threads",
        default=default,
        callback=threads_callback,
        type=click.IntRange(
            1,
        ),
    )


def test_option():
    """Add test option"""

    def test_callback(ctx, param, value):  # pylint: disable=unused-argument
        """Test callback"""
        if value:
            logger.info("Test mode enabled")
            return ["--configfile", "config/config.mimulus.test.yaml"]
        return []

    return click.option(
        "--test",
        help="Run in test mode",
        is_flag=True,
        callback=test_callback,
        expose_value=True,
    )


dry_run_option = click.option("--dry-run", "-n", is_flag=True, help="dry run")
output_file_option = click.option("--output-file", "-o", help="output file name")
