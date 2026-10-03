"""Locate a SurveySimulator (ossssim) checkout used by this project."""
from __future__ import annotations

import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SSIM_NAMES = ("SurveySimulator", "surveysimulator", "Survey-Simulator")


def is_ssim_root(path: Path) -> bool:
    return (path / "src" / "ossssim").is_dir() and (path / "F95").is_dir()


def ssim_root_candidates() -> list[Path]:
    env = os.environ.get("SSIM_ROOT") or os.environ.get("SURVEYSIMULATOR_ROOT")
    cands: list[Path] = []
    if env:
        cands.append(Path(env).expanduser())
    parent = PROJECT_ROOT.parent
    for name in _SSIM_NAMES:
        cands.append(parent / name)
    cands.append(Path("/workspace"))
    cands.extend(Path(__file__).resolve().parents)
    return cands


def find_ssim_root() -> Path:
    """Return the SSim clone that contains ``src/ossssim`` and ``F95``.

    Search order: ``SSIM_ROOT`` / ``SURVEYSIMULATOR_ROOT``, a sibling
    ``SurveySimulator`` directory, ``/workspace``, then parents of this file.
    """
    seen: set[Path] = set()
    for raw in ssim_root_candidates():
        try:
            path = raw.resolve()
        except OSError:
            continue
        if path in seen:
            continue
        seen.add(path)
        if is_ssim_root(path):
            return path
    raise FileNotFoundError(
        "Could not find SurveySimulator (need src/ossssim and F95). "
        "Clone ijiraq/SurveySimulator next to this repo or set SSIM_ROOT."
    )


def add_ssim_src(root: Path | None = None) -> Path:
    """Put ``<ssim>/src`` on ``sys.path`` so ``import ossssim`` works from source."""
    ssim = root or find_ssim_root()
    src = str(ssim / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    return ssim


def default_l7_model() -> Path:
    return find_ssim_root() / "F95" / "tests" / "Models" / "L7model-3.0-9.0"


def load_ossssim_grid_bias():
    """Load ``grid_bias.py`` from SSim without importing the ``ossssim`` package.

    Importing ``ossssim`` pulls ``ossssimlib`` (the Fortran extension). Tests and
    the Sample A wrapper only need the pure-Python helpers.
    """
    import importlib.util

    path = find_ssim_root() / "src" / "ossssim" / "grid_bias.py"
    name = "ossssim_grid_bias_impl"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load ossssim.grid_bias from {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod
