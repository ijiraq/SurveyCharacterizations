#!/usr/bin/env python
"""Run the cold model through pointings.list + MPC.eff and compare with "H < h_lim".

Model rows (OSSOS 1.0 coldl, coldh, kernel; i_free < ifree_max, Huang et al.
2022) get a uniform H in
the survey filter (r) around h_lim and are simulated with ossssim
(comp = cold, colours from the root colour.toml). The detected (flag >= 4) fraction against H should be a
step from 1 to 0 at h_lim: the analytic efficiency survey_debias uses is
f_cut * 1(H < h_lim), where the cut fraction f_cut is applied separately.

Writes MPC_ossssim_validation.{txt,png}.

Usage::

    python validate_ossssim.py --model ../../../KuiperBeltModels/OSSOS/Classical-ModelUsed.dat
"""
from __future__ import annotations

import argparse
import contextlib
import os
import sys
import tomllib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COLD = ("coldl", "coldh", "kernel")


def read_model(path: Path, ifree_max: float):
    """(a, e, i, node, peri, M, epoch) of the cold rows with i_free < ifree_max."""
    from survey_debias.survey_efficiency import ModelPopulation

    m = ModelPopulation.from_file(path, COLD)
    keep = m.ifree < ifree_max
    return np.column_stack([m.a, m.e, m.i, m.node, m.peri, m.M, m.epoch])[keep]


@contextlib.contextmanager
def quiet():
    """Silence the Fortran library's stdout chatter."""
    fd = sys.stdout.fileno()
    saved = os.dup(fd)
    with open(os.devnull, "w") as null:
        os.dup2(null.fileno(), fd)
        try:
            yield
        finally:
            os.dup2(saved, fd)
            os.close(saved)


def main(argv=None) -> int:
    with (HERE / "MPC.toml").open("rb") as fh:
        spec = tomllib.load(fh)
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--model", required=True)
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--ifree-max", type=float, default=float(spec["cut"]["ifree_max"]))
    p.add_argument("--h-lim", type=float, default=float(spec["survey"]["h_lim"]))
    p.add_argument("--half-range", type=float, default=0.75)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args(argv)

    from astropy import units as u
    from ossssim import OSSSSim
    from survey_debias.colours import ColourMap, load_colours

    photspec, _ = load_colours(ROOT / spec["survey"]["colour_file"])
    model_band = spec["survey"]["model_band"]
    colour = ColourMap(photspec, model_band).colour(spec["survey"]["filter"], "cold")
    rng = np.random.default_rng(args.seed)
    model = read_model(Path(args.model), args.ifree_max)
    pick = model[rng.integers(0, len(model), args.n)]
    h_m = rng.uniform(args.h_lim - args.half_range, args.h_lim + args.half_range, args.n)

    sim = OSSSSim(str(HERE), seed=args.seed)
    found = np.zeros(args.n, dtype=bool)
    mu = np.full(args.n, np.nan)
    for k, (a, e, inc, node, peri, M, epoch) in enumerate(pick):
        row = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg, peri=peri * u.deg,
                   M=M * u.deg, H=(h_m[k] - colour) * u.mag, epoch=epoch * u.day, comp="cold")
        with quiet():
            out = sim.simulate(row, colors=photspec, model_band=model_band)
        found[k] = int(out["flag"]) >= 4
        if found[k]:
            mu[k] = float(out["m_int"].value) - h_m[k]
        if k and k % 5000 == 0:
            sys.stderr.write(f"  {k} simulated\n")

    edges = np.arange(args.h_lim - args.half_range, args.h_lim + args.half_range + 1e-9, 0.05)
    centres = 0.5 * (edges[1:] + edges[:-1])
    frac = np.array([found[(h_m >= lo) & (h_m < hi)].mean() for lo, hi in zip(edges[:-1], edges[1:])])
    below = found[h_m < args.h_lim - 0.15].mean()
    above = found[h_m > args.h_lim + 0.15].mean()
    # Edge: H where the detected fraction crosses half of its plateau.
    half = 0.5 * below
    k = int(np.flatnonzero(frac < half)[0])
    h50 = float(np.interp(half, [frac[k], frac[k - 1]], [centres[k], centres[k - 1]]))
    lo10 = float(centres[np.flatnonzero(frac < 0.9 * below)[0]])
    hi90 = float(centres[np.flatnonzero(frac < 0.1 * below)[0]])
    lines = [
        f"# ossssim through pointings.list + MPC.eff: {args.n} cold model rows (i_free < {args.ifree_max:g}),",
        f"# H uniform in {args.h_lim:g} +/- {args.half_range:g}; colour {spec["survey"]["filter"]} - {model_band} = {colour:+.2f}",
        f"detected fraction for H < h_lim - 0.15: {below:.4f}",
        f"detected fraction for H > h_lim + 0.15: {above:.4f}",
        f"50% edge: H = {h50:.3f} (h_lim {args.h_lim:g}; offset {h50 - args.h_lim:+.3f})",
        f"fraction falls below 90% / 10% of the plateau at H = {lo10:.3f} / {hi90:.3f}",
        "",
        "H_centre detected_fraction n",
    ]
    for c, f, lo, hi in zip(centres, frac, edges[:-1], edges[1:]):
        lines.append(f"{c:.3f} {f:.4f} {int(((h_m >= lo) & (h_m < hi)).sum())}")
    (HERE / "MPC_ossssim_validation.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:6]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.step(centres, frac, where="mid", label="ossssim (pointings.list + MPC.eff)")
    ax.plot([edges[0], args.h_lim, args.h_lim, edges[-1]], [1, 1, 0, 0], "k--", lw=0.8,
            label=f"H < {args.h_lim:g}")
    ax.set_xlabel("H (survey filter)")
    ax.set_ylabel("detected fraction (i_free < 4 cold model)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "MPC_ossssim_validation.png", dpi=130)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
