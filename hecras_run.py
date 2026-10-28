#!/usr/bin/env python3
"""Command-line entry point kept for existing commands; the code is in
hecras_linux/run.py (installed as `hecras-run`)."""
import sys

from hecras_linux import run

if __name__ == "__main__":
    sys.exit(run.main())
else:
    sys.modules[__name__] = run
