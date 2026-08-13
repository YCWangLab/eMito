"""Public command-line entry point for eMito."""

from __future__ import annotations

import shutil
import sys
from typing import Callable, Optional, Sequence

from . import __version__, pipeline, reporting, stage_reporting, taxonomy_reporting


HELP = f"""eMito {__version__}

Taxonomy-aware mitochondrial capture probe design.

Usage:
  emito <command> [options]

Commands:
  run               Run the configurable probe-design pipeline
  validate          Validate input files without producing an output directory
  summarize         Count generation, access, collapse, and final probe sets
  stage-summary     Produce detailed probe and k-mer stage statistics
  taxonomy-summary  Count input genomes, species, genera, and families
  info              Show the eMito version and MAFFT availability

Use `emito <command> --help` for command-specific options.
"""


def _dispatch(function: Callable[[Optional[Sequence[str]]], int], args: Sequence[str]) -> int:
    return function(args)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help"}:
        print(HELP)
        return 0
    if args[0] in {"-V", "--version"}:
        print(__version__)
        return 0

    command, command_args = args[0], args[1:]
    if command == "run":
        return _dispatch(pipeline.main, command_args)
    if command == "validate":
        return _dispatch(pipeline.main, [*command_args, "--validate-only"])
    if command == "summarize":
        return _dispatch(reporting.main, command_args)
    if command == "stage-summary":
        return _dispatch(stage_reporting.main, command_args)
    if command == "taxonomy-summary":
        return _dispatch(taxonomy_reporting.main, command_args)
    if command == "info":
        if command_args:
            print("ERROR: `emito info` does not take arguments", file=sys.stderr)
            return 2
        mafft = shutil.which("mafft")
        print(f"eMito: {__version__}")
        print(f"Python: {sys.version.split()[0]}")
        print(f"MAFFT: {mafft or 'not found'}")
        return 0

    print(f"ERROR: unknown command {command!r}\n", file=sys.stderr)
    print(HELP, file=sys.stderr)
    return 2
