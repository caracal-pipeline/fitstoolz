"""The `fitstoolz` command-line apps, one module per subcommand.

Each is a `@shinobi.pystep` plus a `click.Command` built from it by
`_cli.make_command`; `main.py` assembles them into the one lazy group.
"""

from ._outputs import FitsOutputs, outfits_name

__all__ = ["FitsOutputs", "outfits_name"]
