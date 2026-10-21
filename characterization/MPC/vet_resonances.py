#!/usr/bin/env python
"""Resonance vetting of the bright MPC cold-classical candidates (writes MPC.vetting).

Each candidate (main-belt, i_free < ifree_max, q >= q_min, a in the cold range,
H <= h_max in MPC.csv) is integrated with two clones at a -/+ da under the
Sun and the four giant planets (JPL Horizons states at the MPCORB epoch,
WHFast). For every Neptune mean-motion resonance p:q (p - q <= max_order,
p <= max_p) within 0.4 au of the particle's mean a, the resonant angle

    phi = p lambda - q lambda_N - (p - q) varpi

is tested in windows of the integration: it librates in a window if it
avoids some arc of at least ``gap`` degrees. An object is "res" when its
nominal orbit librates in at least half of the windows (as in Gladman et
al. 2008, the clones only flag the classification as insecure).

Usage (from this directory)::

    python vet_resonances.py --t-max 1e7       # ~10 Myr, writes MPC.vetting
    python make_mpc.py                          # applies it to MPC.csv
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

from make_mpc import COLD_A_MIN, INNER_EDGE, OUTER_EDGE, read_csv

HERE = Path(__file__).resolve().parent
PLANETS = ("Jupiter", "Saturn", "Uranus", "Neptune")
A_NEPTUNE = 30.07
INNER_MASS = 6.0e-6        # Mercury to Mars plus the Moon, added to the Sun [Msun]


def candidates(rows, args):
    out = []
    for r in rows:
        if not (INNER_EDGE <= r["a"] <= OUTER_EDGE and args.a_min <= r["a"] <= args.a_max):
            continue
        if not r["ifree"] < args.ifree_max or r["q"] < args.q_min or r["h_mpc"] > args.h_max:
            continue
        out.append(r)
    return out


def resonances(a_mean: float, max_order: int, max_p: int, width: float = 0.4):
    out = []
    for p in range(2, max_p + 1):
        for k in range(1, max_order + 1):
            q = p - k
            if q < 1 or math.gcd(p, q) != 1:
                continue
            a_res = A_NEPTUNE * (p / q) ** (2.0 / 3.0)
            if abs(a_res - a_mean) <= width:
                out.append((p, q, a_res))
    return out


def librating(phi_deg: np.ndarray, gap: float) -> bool:
    """True if phi avoids an arc of at least ``gap`` degrees."""
    s = np.sort(np.mod(phi_deg, 360.0))
    gaps = np.diff(np.concatenate([s, [s[0] + 360.0]]))
    return bool(gaps.max() >= gap)


def build(rows, epoch_jd: float, da: float, mpc_epoch_jd):
    import rebound

    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    date = f"JD{epoch_jd:.1f}"
    sim.add("Sun", date=date)
    sim.particles[0].m += INNER_MASS
    for name in PLANETS:
        sim.add(name, date=date)
    sun = sim.particles[0]
    gm = sim.G * sun.m
    index = []
    for r in rows:
        dt_yr = (epoch_jd - mpc_epoch_jd(r["epoch"])) / 365.25
        for k, shift in enumerate((0.0, -da, +da)):
            a = r["a"] + shift
            n = math.sqrt(gm / a ** 3)
            mean_anom = math.radians(float(r["M"])) + n * dt_yr
            sim.add(primary=sun, m=0.0, a=a, e=r["e"], inc=math.radians(r["i"]),
                    Omega=math.radians(float(r["node"])), omega=math.radians(float(r["peri"])),
                    M=mean_anom)
            index.append((r["desig"], k))
    sim.N_active = 1 + len(PLANETS)
    sim.testparticle_type = 0
    sim.move_to_com()
    sim.integrator = "whfast"
    sim.dt = 0.5
    return sim, index


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--csv", default=str(HERE / "MPC.csv"))
    p.add_argument("--out", default=str(HERE / "MPC.vetting"))
    p.add_argument("--h-max", type=float, default=6.0)
    p.add_argument("--ifree-max", type=float, default=4.0)
    p.add_argument("--q-min", type=float, default=37.0)
    p.add_argument("--a-min", type=float, default=COLD_A_MIN)
    p.add_argument("--a-max", type=float, default=OUTER_EDGE)
    p.add_argument("--t-max", type=float, default=1.0e7, help="years")
    p.add_argument("--n-out", type=int, default=20000)
    p.add_argument("--windows", type=int, default=10)
    p.add_argument("--gap", type=float, default=20.0, help="degrees")
    p.add_argument("--da", type=float, default=0.005, help="clone offset in a [au]")
    p.add_argument("--max-order", type=int, default=6)
    p.add_argument("--max-p", type=int, default=40)
    args = p.parse_args(argv)

    from astropy.time import Time

    rows, _ = read_csv(Path(args.csv))
    cand = candidates(rows, args)
    epochs = sorted({r["epoch"] for r in cand})
    epoch_jd = Time(max(epochs, key=lambda e: sum(r["epoch"] == e for r in cand))).jd

    def mpc_epoch_jd(text):
        return Time(text).jd

    sim, index = build(cand, epoch_jd, args.da, mpc_epoch_jd)
    n_p = len(index)
    sys.stderr.write(f"{len(cand)} candidates ({n_p} particles), epoch JD {epoch_jd}, "
                     f"{args.t_max:.3g} yr\n")
    times = np.linspace(0.0, args.t_max, args.n_out)
    a = np.zeros((args.n_out, n_p))
    lam = np.zeros((args.n_out, n_p))
    pom = np.zeros((args.n_out, n_p))
    lam_n = np.zeros(args.n_out)
    t0 = time.time()
    first = 1 + len(PLANETS)
    for j, t in enumerate(times):
        sim.integrate(t, exact_finish_time=0)
        sun = sim.particles[0]
        nep = sim.particles[len(PLANETS)].orbit(primary=sun)
        lam_n[j] = nep.l
        for k in range(n_p):
            o = sim.particles[first + k].orbit(primary=sun)
            a[j, k], lam[j, k], pom[j, k] = o.a, o.l, o.pomega
        if j and j % (args.n_out // 10) == 0:
            sys.stderr.write(f"  t = {t:.3g} yr ({time.time() - t0:.0f} s)\n")

    windows = np.array_split(np.arange(args.n_out), args.windows)
    results = {}
    for k, (desig, clone) in enumerate(index):
        a_mean = float(np.mean(a[:, k]))
        if not np.all(np.isfinite(a[:, k])) or a[:, k].max() > 1000 or a[:, k].min() < 0:
            results.setdefault(desig, []).append((clone, a_mean, None, 0.0))
            continue
        best = (None, 0.0)
        for pp, qq, a_res in resonances(a_mean, args.max_order, args.max_p):
            phi = np.degrees(pp * lam[:, k] - qq * lam_n - (pp - qq) * pom[:, k])
            frac = np.mean([librating(phi[w], args.gap) for w in windows])
            if frac > best[1]:
                best = (f"{pp}:{qq}", frac)
        results.setdefault(desig, []).append((clone, a_mean, best[0], best[1]))

    lines = [
        "# Resonance vetting of the MPC cold candidates (vet_resonances.py)",
        f"# {len(cand)} objects with H <= {args.h_max:g}, {args.a_min:g} <= a <= {args.a_max:g}, "
        f"i_free < {args.ifree_max:g}, q >= {args.q_min:g}; nominal and a -/+ {args.da:g} au clones;",
        f"# {args.t_max:.3g} yr WHFast, giant planets from JPL Horizons at JD {epoch_jd}; "
        f"libration = phi avoids a {args.gap:g} deg arc,",
        f"# counted in {args.windows} windows. cl = res when the nominal orbit librates in >= half.",
        "# desig cl  H  <a>  nominal_resonance  libration_fraction (nominal, a-da, a+da)  "
        "secure (all three agree)  resonances any orbit librates in",
    ]
    for r in sorted(cand, key=lambda r: r["h"]):
        res = sorted(results[r["desig"]])
        nominal = res[0]
        cl = "res" if nominal[2] and nominal[3] >= 0.5 else "cla"
        fr = " ".join(f"{x[3]:.1f}" for x in res)
        clones = {x[2] for x in res if x[3] >= 0.5}
        agree = all(x[3] >= 0.5 for x in res) or all(x[3] < 0.5 for x in res)
        secure = "secure" if agree else "insecure"
        lines.append(f"{r['desig']:12s} {cl:4s} H={r['h']:.2f} a={nominal[1]:.2f} "
                     f"{nominal[2] or '-'} {fr} {secure} {','.join(sorted(c for c in clones if c))}")
    Path(args.out).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
