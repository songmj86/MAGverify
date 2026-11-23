"""Installed command-line entry point for MAGverify."""
import sys

if sys.version_info < (3, 10):
    sys.exit('MAGverify needs Python 3.10 or newer (found %d.%d)' % sys.version_info[:2])

from scripts.core_magverify_cli import main

__all__ = ["main"]
