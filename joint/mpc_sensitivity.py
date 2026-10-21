#!/usr/bin/env python
"""Sensitivity of the joint cold fit with the MPC complete survey (cold_all_mpc.toml).

Maximum-likelihood fits (tapered, no MCMC) of variants of cold_all_mpc.toml:
the MPC H_r limit, the H_r scatter sigma_h, the H_r zero point (+/- 0.1),
and the vetting choice for 2025 ME398 (one-opposition, a = 41.7 au; counted
as cold in that variant). Writes joint_mpc_sensitivity/summary.{json,txt}.

    cd joint && python mpc_sensitivity.py [--jobs 8]
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import tomllib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CONFIG = HERE / "cold_all_mpc.toml"
OUT = HERE / "joint_mpc_sensitivity"


def variants(cfg: dict, root: Path) -> dict:
    out = {"base": cfg}
    for h_lim in (5.0, 5.25):
        v = copy.deepcopy(cfg)
        v["survey"][0]["h_lim"] = h_lim
        out[f"h_lim={h_lim:g}"] = v
    for sigma_h in (0.0, 0.15):
        v = copy.deepcopy(cfg)
        v["survey"][0]["sigma_h"] = sigma_h
        out[f"sigma_h={sigma_h:g}"] = v
    src = root / "characterization" / "MPC" / "MPC.csv"
    lines = src.read_text().splitlines()
    head = [x for x in lines if x.startswith("#")]
    rows = list(csv.DictReader([x for x in lines if not x.startswith("#")]))

    def variant_csv(name, edit):
        new = [edit(dict(r)) for r in rows]
        path = OUT / f"MPC_{name}.csv"
        with path.open("w", newline="") as fh:
            fh.write("\n".join(head) + "\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(new)
        v = copy.deepcopy(cfg)
        v["survey"][0]["detections"] = str(path)
        return v

    # H_r zero point: OSSOS/CFEPS Hsur runs 0.11 mag brighter (MPC_calibration.txt).
    for shift in (-0.1, 0.1):
        out[f"H_r {shift:+g}"] = variant_csv(
            f"shift{shift:+g}", lambda r, s=shift: {**r, "h": f"{float(r['h']) + s:.3f}"})

    def with_me398(r):
        if r["desig"] == "2025ME398":
            r["cl"], r["p"], r["comp"] = "cla", "m", "cold"
        return r

    out["with 2025 ME398"] = variant_csv("with_2025ME398", with_me398)
    return out


def run_one(item):
    label, cfg = item
    from survey_debias import hfit
    from survey_debias.joint import fit_joint
    from survey_debias.joint_cli import build_terms, summarize_counts

    like, _ = build_terms(cfg, HERE)
    fit = fit_joint(like)
    shape, h_o, _ = like.unpack(fit.theta)
    n_lt = {f"N(H_r<{h:g})": float(hfit.cumulative(h, *shape, h_o)) for h in (4.5, 5.0, 5.5, 6.0)}
    rel = fit.profile_loglike - np.nanmax(fit.profile_loglike)
    inside = fit.profile_alpha[rel >= -0.5]
    mpc = like.surveys[0]
    return label, {
        "mle": dict(zip(like.names, map(float, fit.theta))), "loglike": float(fit.loglike),
        "counts": summarize_counts(like, fit.theta), **n_lt,
        "mpc": {"n": mpc.n, "h_lim": mpc.h_lim, "sigma_h": mpc.sigma_h, "colour": mpc.colour},
        "alpha_within_0p5": [float(inside.min()), float(inside.max())],
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--jobs", type=int, default=8)
    args = p.parse_args(argv)
    OUT.mkdir(exist_ok=True)
    with CONFIG.open("rb") as fh:
        cfg = tomllib.load(fh)
    root = (HERE / cfg["joint"]["root"]).resolve()
    todo = variants(cfg, root)
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        results = dict(pool.map(run_one, todo.items()))
    results = {k: results[k] for k in todo}
    (OUT / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    lines = ["MPC sensitivity: maximum likelihood (tapered) of cold_all_mpc.toml variants", "",
             f"{'variant':18s} {'n_MPC':>5s} {'exp_MPC':>7s} {'alpha':>6s} {'beta':>6s} "
             f"{'H_B':>6s} {'N(<5)':>6s} {'N(<5.5)':>7s} {'N(<6)':>6s} {'lnL':>8s}"]
    for label, r in results.items():
        m = r["mle"]
        lines.append(f"{label:18s} {r['mpc']['n']:5d} {r['counts']['MPC']['expected']:7.2f} "
                     f"{m['alpha']:6.3f} {m['beta']:6.3f} {m['h_b']:6.2f} "
                     f"{r['N(H_r<5)']:6.2f} {r['N(H_r<5.5)']:7.2f} {r['N(H_r<6)']:6.2f} "
                     f"{r['loglike']:8.2f}")
    lines += ["", "N(H_r<H): the whole cold population (all inclinations) at the MLE;",
              "exp_MPC: expected MPC cold count (i < 5 deg, below the limit) at the MLE."]
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
