"""JWST Sample A grid-cell helpers.

Loads SurveySimulator ``src/ossssim/grid_bias.py`` without importing the
``ossssim`` package (that would pull ossssimlib). Module-level field defaults
are this survey. Set ``SSIM_ROOT`` or clone SurveySimulator next to this repo.
"""
from __future__ import annotations

import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from ssim_paths import load_ossssim_grid_bias  # noqa: E402

_mod = load_ossssim_grid_bias()
globals().update({k: v for k, v in vars(_mod).items() if not k.startswith("__")})
SURVEY = JWST_SAMPLE_A  # noqa: F405
