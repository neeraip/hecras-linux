#!/usr/bin/env python3
"""Command-line entry point kept for existing commands; the code is in
hecras_linux/preprocess.py (installed as `hecras-preprocess`)."""
import sys

from hecras_linux import preprocess

if __name__ == "__main__":
    preprocess.main()
else:
    sys.modules[__name__] = preprocess
