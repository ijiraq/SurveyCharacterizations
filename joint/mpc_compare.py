#!/usr/bin/env python
"""Joint fits with and without the MPC complete survey, and MPC posterior-predictive checks.

For each posterior (joint_samples.csv of survey_debias_joint, or the
hfit_samples.csv of a single-survey survey_debias_hfit run: columns alpha,
beta, h_b, h_o in the model band) this reports N(H_r < H) at the bright end
and the number of cold MPC objects it predicts below the MPC limit,

    Lambda_MPC = f_cut int lambda_M(H) 1(H < h_lim) dH

(lambda_M: the LF shifted by the MPC colour and spread by sigma_c and
sigma_h, as in the joint fit), against the observed count: the fraction of
posterior-predictive draws N ~ Poisson(Lambda_MPC) at or below, and at or
above, the observed n.

    cd joint && python mpc_compare.py [--hfit-samples path/hfit_samples.csv ...]
    python mpc_compare.py --config coldfree_all_mpc.toml --without joint_coldfree_all \
        --with joint_coldfree_all_mpc --out joint_mpc_compare_coldfree

Writes joint_mpc_compare/compare.{txt,png}.
"""
from __future__ import annotations

import argparse
import copy
import json
import tomllib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = {"without MPC": HERE / "joint_cold_all_nompc", "with MPC": HERE / "joint_cold_all_mpc"}
OUT = HERE / "joint_mpc_compare"
H_POINTS = (4.5, 5.0, 5.5, 6.0)


def mpc_only(config: Path):
    """JointLikelihood holding just the MPC term of a configuration."""
    from survey_debias.joint_cli import build_terms

    with config.open("rb") as fh:
        cfg = tomllib.load(fh)
    cfg = copy.deepcopy(cfg)
    cfg["survey"] = [s for s in cfg["survey"] if s.get("format") == "complete"]
    cfg["joint"]["priority"] = [s["name"] for s in cfg["survey"]]
    like, _ = build_terms(cfg, config.parent)
    return like


def read_samples(path: Path) -> tuple[list[str], np.ndarray]:
    names = path.read_text().splitlines()[0].split(",")
    return names, np.loadtxt(path, delimiter=",", skiprows=1)


def bright_counts(samples: np.ndarray, names: list[str]) -> dict:
    from survey_debias import hfit

    col = {n: k for k, n in enumerate(names)}
    a, b, hb, ho = (samples[:, col[n]] for n in ("alpha", "beta", "h_b", "h_o"))
    return {h: hfit.cumulative(h, a, b, hb, ho) for h in H_POINTS}


def predicted_mpc(like, samples: np.ndarray, names: list[str], n_draw: int = 4000,
                  seed: int = 0) -> np.ndarray:
    col = {n: k for k, n in enumerate(names)}
    term = like.surveys[0]
    key = f"delta_{term.filter}"
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(samples), min(n_draw, len(samples)), replace=False)
    lam = np.empty(len(pick))
    for j, k in enumerate(pick):
        s = samples[k]
        shape = (s[col["alpha"]], s[col["beta"]], s[col["h_b"]])
        delta = {term.filter: s[col[key]] if key in col else 0.0}
        lam[j] = like.expected(term, shape, s[col["h_o"]], delta)[0]
    return lam


def q(x) -> str:
    lo, med, hi = np.percentile(x, [16, 50, 84])
    return f"{med:.3g} (-{med - lo:.2g} +{hi - med:.2g})"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--config", default=str(HERE / "cold_all_mpc.toml"))
    p.add_argument("--without", default=str(RUNS["without MPC"]),
                   help="run directory of the fit without the MPC")
    p.add_argument("--with", dest="with_mpc", default=str(RUNS["with MPC"]),
                   help="run directory of the fit with the MPC")
    p.add_argument("--out", default=str(OUT), help="output directory")
    p.add_argument("--hfit-samples", action="append", default=[],
                   help="hfit_samples.csv of a survey_debias_hfit run (repeatable)")
    args = p.parse_args(argv)

    like = mpc_only(Path(args.config))
    term = like.surveys[0]
    n_obs = term.n
    f_cut = float(term.eps.meta["f_cut"])
    runs = {"without MPC": Path(args.without), "with MPC": Path(args.with_mpc)}
    posteriors = {label: d / "joint_samples.csv" for label, d in runs.items()
                  if (d / "joint_samples.csv").is_file()}
    out_dir = Path(args.out)
    for path in args.hfit_samples:
        posteriors[f"hfit {Path(path).parent.name}"] = Path(path)

    rng = np.random.default_rng(1)
    lines = [f"MPC: {n_obs} cold objects with {term.filter} < {term.h_lim:g} "
             f"(f_cut {f_cut:.3f}, colour {term.colour:+.2f}, sigma_h {term.sigma_h:g})", ""]
    curves = {}
    for label, path in posteriors.items():
        names, samples = read_samples(path)
        counts = bright_counts(samples, names)
        lam = predicted_mpc(like, samples, names)
        draws = rng.poisson(lam)
        lines.append(f"{label} ({path.relative_to(HERE) if path.is_relative_to(HERE) else path}):")
        summary = path.parent / "joint_summary.json"
        if summary.is_file():
            post = json.loads(summary.read_text()).get("posterior", {})
            lines.append("  " + ", ".join(f"{k} {v['median']:.3f} -{v['minus']:.3f} "
                                          f"+{v['plus']:.3f}" for k, v in post.items()
                                          if k in ("alpha", "beta", "h_b")))
        for h, x in counts.items():
            lines.append(f"  N(H_r < {h:g}) = {q(x)}")
        lines.append(f"  predicted MPC count {q(lam)}; P(N <= {n_obs}) = "
                     f"{np.mean(draws <= n_obs):.3f}, P(N >= {n_obs}) = "
                     f"{np.mean(draws >= n_obs):.3f}")
        lines.append("")
        curves[label] = (names, samples)
    text = "\n".join(lines) + "\n"
    out_dir.mkdir(exist_ok=True)
    (out_dir / "compare.txt").write_text(text)
    print(text)
    plot(curves, term, f_cut, out_dir)
    return 0


def plot(curves, term, f_cut, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from survey_debias import hfit

    hg = np.linspace(3.5, 8.0, 200)
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    colours = {"without MPC": "#c45c26", "with MPC": "#1f4e79"}
    for label, (names, samples) in curves.items():
        col = {n: k for k, n in enumerate(names)}
        pick = samples[np.random.default_rng(0).choice(len(samples), min(400, len(samples)),
                                                       replace=False)]
        cum = np.array([hfit.cumulative(hg, s[col["alpha"]], s[col["beta"]], s[col["h_b"]],
                                        s[col["h_o"]]) for s in pick])
        lo, med, hi = np.percentile(cum, [16, 50, 84], axis=0)
        c = colours.get(label, "0.4")
        ax.fill_between(hg, lo, hi, color=c, alpha=0.25)
        ax.plot(hg, med, color=c, lw=1.8, label=f"{label} (median, 68%)")
    h_r = np.sort(term.like.h - term.colour)
    ax.step(h_r, np.arange(1, len(h_r) + 1) / f_cut, where="post", color="k", lw=1.2,
            marker="o", ms=3, label=f"MPC: {term.filter} < {term.h_lim:g}, / f_cut")
    ax.axvline(term.h_lim - term.colour, color="k", ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(0.1, 3e3)
    ax.set_xlabel(r"$H_r$")
    ax.set_ylabel(r"$N(<H)$, whole cold population")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / "compare.png", dpi=150)


if __name__ == "__main__":
    raise SystemExit(main())
