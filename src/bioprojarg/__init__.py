"""
Bioproject ARG inference workflow
"""

__version__ = "undefined"

try:
    from . import _version  # pylint: disable=import-self

    __version__ = _version.version
except ImportError:
    pass
