"""Compatibility module for the nanobind LillyMol binding.

The nanobind changeover exposes these symbols from the main `lillymol` module.
This module keeps older split-module imports working.
"""

from lillymol import *  # noqa: F401,F403
