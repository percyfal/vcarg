"""Environment for BioprojArg"""

from . import __version__


class Environment:  # pylint: disable=too-few-public-methods
    """Environment for BioprojArg"""

    def __init__(self):
        self.debug = False
        self.version = __version__
