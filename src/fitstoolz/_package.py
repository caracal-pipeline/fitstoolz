"""Package-level constants, kept out of the package initializer.

An `__init__.py` that holds implementation is what RUF067
(`non-empty-init-module`) is about: importing the package then runs it, and
`from fitstoolz import X` says nothing about where `X` is defined.
"""

import logging
from importlib import metadata

__version__ = metadata.version(__package__)

#: The logger every module in the package emits through, directly or (as
#: `reader.py` does) through a `__name__` child of it.
LOGGER = "fitstoolz"

# Library convention, the same one shinobi keeps for `shinobi.*`: modules
# only ever *emit* through this logger and never attach a handler, so
# importing fitstoolz configures nothing. The NullHandler is what makes that
# silence deliberate rather than leaving logging's last-resort stderr echo of
# unhandled WARNING+ records to decide. The console handler belongs to the
# application -- see `fitstoolz.apps._cli`.
logging.getLogger(LOGGER).addHandler(logging.NullHandler())
