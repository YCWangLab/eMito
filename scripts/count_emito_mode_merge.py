#!/usr/bin/env python3
"""Compatibility wrapper for ``emito summarize``."""

from emito.reporting import main


if __name__ == "__main__":
    raise SystemExit(main())
