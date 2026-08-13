#!/usr/bin/env python3
"""Compatibility wrapper for ``emito stage-summary``."""

from emito.stage_reporting import main


if __name__ == "__main__":
    raise SystemExit(main())
