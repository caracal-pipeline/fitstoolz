"""fitstoolz's test suite.

`InitTest` builds the example FITS files the app tests run against and
cleans up what it made; `TESTDIR` is where it puts them.
"""

from ._helpers import TESTDIR, InitTest

__all__ = ["TESTDIR", "InitTest"]
