"""Build a `click.Command` from a `@shinobi.pystep` StepRef.

Replaces scabha's `clickify_parameters` + YAML parser configs: options are
derived from the step's pydantic ``inputs_model`` by shinobi's
``build_options`` (dtype, choices, abbreviations, bool/list handling), so a
step's signature is the single schema authority. Mirrors simms 3.0's
``simms.apps.main._make_command``.

This module is also fitstoolz's *application* layer, and so the one place
that configures logging. The package itself only ever calls
``logging.getLogger`` (see `fitstoolz._package`), the same library/application
split shinobi keeps -- ``shinobi.logsetup`` attaches shinobi's file handler
from the CLI and nothing in the package touches a handler. ``set_logger``
used to live in fitstoolz's package initializer and did both jobs at once.
"""

from __future__ import annotations

import logging

import click
from shinobi.clickutil import build_options, unflatten_kwargs
from shinobi.steps.dispatch import _dispatch

from fitstoolz import LOGGER


def configure_logging(level: str | int = "INFO") -> None:
    """Attach fitstoolz's console handler, at `level`.

    Called once per command invocation, from the callback below, so every
    subcommand is configured -- `slice`, `unstack` and `sanitise` declare
    ``--log-level`` like the rest but never called ``set_logger``, so the
    level they were given did nothing.

    An unrecognised name lands on INFO. It used to land on 10 (DEBUG), so a
    typo turned the logging all the way *up* -- the same fallback simms
    corrected in shinobi-dosho/simms#197.

    Args:
        level: A level name (``"INFO"``, ``"DEBUG"``, ...), or a logging level.
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    logger = logging.getLogger(LOGGER)
    logger.setLevel(level)

    if not any(not isinstance(handler, logging.NullHandler) for handler in logger.handlers):
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s-%(name)s-%(levelname)-8s| %(message)s"))
        logger.addHandler(handler)
    for handler in logger.handlers:
        handler.setLevel(level)


def make_command(step, *, positional: str | None = None, extra_options=()) -> click.Command:
    """Build a `click.Command` for a `@shinobi.pystep` StepRef.

    Args:
        step: The `StepRef` produced by `@shinobi.pystep`.
        positional: Name of the input field to render as a `click.Argument`
            rather than an option (``build_options`` only emits ``--options``).
            This is the scabha ``policies.positional`` equivalent.
        extra_options: Extra `click.Parameter`s to prepend (eager flags etc.).

    Returns:
        A `click.Command` whose callback re-nests the flat kwargs and
        dispatches the step in-process via shinobi, exactly as
        `shinobi.cli`'s ``run`` command does. ``log_level`` is dropped from
        the options and taken from the root group instead.
    """
    model = step.step.inputs_model
    options = [opt for opt in build_options(model) if opt.name != "log_level"]

    params = list(extra_options)
    for opt in options:
        if opt.name == positional:
            params.append(click.Argument([positional], required=True, type=opt.type))
        else:
            params.append(opt)

    def _callback(**raw):
        ctx = click.get_current_context()
        kwargs = unflatten_kwargs(model, raw)
        kwargs["log_level"] = ctx.obj["log_level"]
        # Before dispatch, so anything the step logs on its way in is carried.
        configure_logging(kwargs["log_level"])
        result = _dispatch(step.step, step.func, **kwargs)
        if not result.success:
            raise click.ClickException(f"{step.step.name!r} failed (returncode {result.returncode}).")

    return click.Command(
        name=step.step.name,
        params=params,
        callback=_callback,
        help=step.step.info,
        no_args_is_help=True,
    )
