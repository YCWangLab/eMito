#!/usr/bin/env python3
"""Compatibility wrapper for ``emito run``."""

from emito.pipeline import main


if __name__ == "__main__":
    raise SystemExit(main())
