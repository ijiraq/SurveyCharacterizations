#!/usr/bin/env python
"""Write pointings.list and MPC_<lon>.eff: an ossssim imitation of "complete for H < h_lim".

OSSSSim efficiencies depend on apparent magnitude and rate, never on H. At
opposition the rate of a TNO fixes its distance r, so a magnitude limit
that depends on rate detects an object when its H (in the survey filter,
MPC.toml: r) is below h_lim:

    rate(r)  = C (v_E - v_1 / sqrt(r)) / (r - R_E)          (circular orbit)
    m_lim(r) = h_lim + 5 log10(r (r - R_E)) + phase

with R_E and v_E Earth's heliocentric distance and tangential speed on the
tile's date, v_1 = 29.785 km/s, C = 4.9637 "/hr per (km/s / au), and phase
the H-G (G = 0.15) term ossssim adds at the small phase angles inside a tile
(about 0.03 mag; averaged over the tile at the cold belt's latitudes).

pointings.list tiles the ecliptic band |beta| <= beta_max with polygons
dlon wide in ecliptic longitude, each observed at the date its centre is at
opposition, from the geocentre (500), fill 1, with its own MPC_<lon>.eff.
Neighbouring tiles overlap by ``margin`` deg so an object that moves
between two tiles' dates is not missed; efficiency is 0 or 1, so seeing an
object twice changes nothing (never run survey_debias_eff on these tiles:
it needs disjoint tiles). A tile that crosses RA = 0 is written twice, at
RA0 and RA0 +/- 360, since the in-field test is planar in RA.

Each .eff has up to 10 rate ranges, one per distance interval. In each, the
efficiency is a narrow single tanh (width ``width`` mag) whose centre
m0 = M0 + alpha_rate * rate is linear in rate and exact at the interval
ends. There is no photometric noise (mag_error 0.001 mag, constant: an
exactly zero first coefficient gives NaN measured magnitudes fainter than
m = 21 in the current ossssim build) and no tracking
loss, so ossssim detects an object when its true magnitude is below
m0(rate). MPC H scatter (MPC.toml sigma_h) is not included: add it to the
model's H if it is wanted.

Usage (from this directory)::

    python make_pointings.py               # values from MPC.toml
"""
from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
C_RATE = 206264.806 * 3600.0 / 1.495978707e8    # "/hr per (km/s / au)
V_1AU = 29.7847                                  # circular speed at 1 au [km/s]
R_EDGES = (30.0, 37.0, 39.5, 42.0, 44.5, 47.0, 49.5, 52.0, 56.0, 62.0, 100.0)
G_SLOPE = 0.15


def opposition_rate(r, r_e: float = 1.0, v_e: float = V_1AU):
    r = np.asarray(r, dtype=float)
    return C_RATE * (v_e - V_1AU / np.sqrt(r)) / (r - r_e)


def phase_mag(alpha_rad, g: float = G_SLOPE):
    """ossssim's H-G term: mag - H - 5 log10(r Delta)."""
    t = np.tan(np.asarray(alpha_rad) / 2.0)
    phi1 = np.exp(-3.33 * t ** 0.63)
    phi2 = np.exp(-1.87 * t ** 1.22)
    return -2.5 * np.log10((1.0 - g) * phi1 + g * phi2)


def mean_phase(r: float, half_lon: float, beta_sd: float, r_e: float = 1.0) -> float:
    """H-G term averaged over a tile: uniform in longitude, normal in latitude [deg]."""
    rng = np.random.default_rng(0)
    dlon = np.radians(rng.uniform(-half_lon, half_lon, 20000))
    beta = np.radians(rng.normal(0.0, beta_sd, 20000))
    theta = np.arccos(np.cos(dlon) * np.cos(beta))      # angle from the anti-solar point
    alpha = np.arcsin(np.clip(r_e * np.sin(theta) / r, -1.0, 1.0))
    return float(np.mean(phase_mag(alpha)))


def m_lim(r, h_lim: float, r_e: float = 1.0, phase: float = 0.0):
    r = np.asarray(r, dtype=float)
    return h_lim + 5.0 * np.log10(r * (r - r_e)) + phase


def rate_ranges(h_lim: float, r_e: float = 1.0, v_e: float = V_1AU, phase: float = 0.0,
                edges=R_EDGES):
    """[(rate_lo, rate_hi, M0, alpha, max_err)] from slow to fast; rates in "/hr."""
    out = []
    for r_near, r_far in zip(edges[:-1], edges[1:]):
        v_far = float(opposition_rate(r_far, r_e, v_e))
        v_near = float(opposition_rate(r_near, r_e, v_e))
        m_far = float(m_lim(r_far, h_lim, r_e, phase))
        m_near = float(m_lim(r_near, h_lim, r_e, phase))
        alpha = (m_near - m_far) / (v_near - v_far)
        m0 = m_far - alpha * v_far
        r = np.linspace(r_near, r_far, 201)
        err = float(np.max(np.abs(m0 + alpha * opposition_rate(r, r_e, v_e)
                                  - m_lim(r, h_lim, r_e, phase))))
        out.append((v_far, v_near, m0, alpha, err, r_near, r_far))
    out.sort()
    # Open-ended: slowest range down to 0, fastest up to 20 "/hr.
    out[0] = (0.0,) + out[0][1:]
    out[-1] = out[-1][:1] + (20.0,) + out[-1][2:]
    return out


def write_eff(path: Path, h_lim: float, width: float, filt: str, r_e: float, v_e: float,
              phase: float, jd: float) -> list:
    ranges = rate_ranges(h_lim, r_e, v_e, phase)
    lines = [
        f"# MPC complete survey: detected when H_{filt} is below",
        f"# h_lim = {h_lim:g} for an object at opposition (make_pointings.py).",
        f"# earth: jd= {jd:.5f} r_e= {r_e:.6f} v_e= {v_e:.4f} phase= {phase:.4f}",
        "# Per rate range: single tanh, m0 = M0 + alpha_rate * rate [\"/hr], exact",
        "# at the range's distance edges; max |m0 - m_lim| inside is given.",
        "rate_cut= 0.0 20.0 0.0 180.0",
        "track_frac= 1.0 99.0 0.0",
        "mag_error= 0.001 0.0 99.0 0.0 99.0 0.0",
        "phot_frac= 1.0 0.0 0.0",
        f"filter= {filt}",
    ]
    for lo, hi, m0, alpha, err, r_near, r_far in ranges:
        lines += [
            f"# r from {r_near:g} to {r_far:g} au; max error {err:.3f} mag",
            f"rates= {lo:.4f} {hi:.4f}",
            "function= single",
            f"single_param= 1.0 {m0:.4f} {width:g} {alpha:.5f}",
        ]
    path.write_text("\n".join(lines) + "\n")
    return ranges


def opposition_jd(lon_deg: float, year_start: str) -> float:
    """JD in the year from ``year_start`` when the Sun is at ecliptic longitude lon + 180."""
    from astropy import units as u
    from astropy.coordinates import GeocentricMeanEcliptic, get_sun
    from astropy.time import Time

    t0 = Time(year_start)
    t = t0 + np.arange(0.0, 366.0, 0.25) * u.day
    sun = get_sun(t).transform_to(GeocentricMeanEcliptic(equinox="J2000", obstime=t))
    target = (lon_deg + 180.0) % 360.0
    d = (sun.lon.deg - target + 180.0) % 360.0 - 180.0
    k = int(np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0))[0])
    frac = -d[k] / (d[k + 1] - d[k])
    return float(t[k].jd + frac * 0.25)


def earth_state(jd: float) -> tuple[float, float]:
    """Earth's heliocentric distance [au] and speed perpendicular to the Sun line [km/s]."""
    from astropy import units as u
    from astropy.coordinates import get_body_barycentric_posvel
    from astropy.time import Time

    t = Time(jd, format="jd")
    pe, ve = get_body_barycentric_posvel("earth", t)
    ps, vs = get_body_barycentric_posvel("sun", t)
    r = (pe - ps).xyz.to(u.au).value
    v = (ve - vs).xyz.to(u.km / u.s).value
    rhat = r / np.linalg.norm(r)
    v_t = v - np.dot(v, rhat) * rhat
    return float(np.linalg.norm(r)), float(np.linalg.norm(v_t))


def tile_vertices(lon0: float, lon1: float, beta: float, step: float = 2.0):
    """Polygon (lon, lat) in ecliptic J2000 around [lon0, lon1] x [-beta, beta]."""
    n = max(1, int(math.ceil((lon1 - lon0) / step)))
    lons = np.linspace(lon0, lon1, n + 1)
    bottom = [(lo, -beta) for lo in lons]
    top = [(lo, beta) for lo in lons[::-1]]
    return bottom + top


def eff_name(lon_c: float) -> str:
    return f"MPC_{int(round(lon_c)):03d}.eff"


def write_survey(h_lim: float, filt: str, beta: float, dlon: float, margin: float,
                 width: float, year_start: str, beta_sd: float) -> tuple[int, float]:
    from astropy import units as u
    from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord

    for old in HERE.glob("MPC*.eff"):
        old.unlink()
    lines = [
        "# MPC complete survey (make_pointings.py): ecliptic band |beta| <= "
        f"{beta:g} deg in {dlon:g} deg tiles (+/- {margin:g} deg overlap),",
        "# each at its own opposition date with its own .eff; fill 1, geocentre.",
        "# Tiles crossing RA = 0 are written twice (RA0 and RA0 +/- 360).",
        "# poly n RA0 Dec0 JD fill observer eff, then n lines of dx dy [deg],",
        "# dx = (RA - RA0) cos(Dec).",
    ]
    n_tiles, worst = 0, 0.0
    phase = mean_phase(44.0, 0.5 * dlon + margin, beta_sd)
    for lon_c in np.arange(0.5 * dlon, 360.0, dlon):
        verts = tile_vertices(lon_c - 0.5 * dlon - margin, lon_c + 0.5 * dlon + margin, beta)
        ecl = SkyCoord([v[0] for v in verts] * u.deg, [v[1] for v in verts] * u.deg,
                       frame=BarycentricMeanEcliptic(equinox="J2000"))
        eq = ecl.icrs
        centre = SkyCoord(lon_c * u.deg, 0.0 * u.deg,
                          frame=BarycentricMeanEcliptic(equinox="J2000")).icrs
        ra0, dec0 = float(centre.ra.deg), float(centre.dec.deg)
        dra = (eq.ra.deg - ra0 + 180.0) % 360.0 - 180.0
        dx = dra * np.cos(np.radians(eq.dec.deg))
        dy = eq.dec.deg - dec0
        jd = opposition_jd(lon_c, year_start)
        r_e, v_e = earth_state(jd)
        name = eff_name(lon_c)
        ranges = write_eff(HERE / name, h_lim, width, filt, r_e, v_e, phase, jd)
        worst = max(worst, max(r[4] for r in ranges[1:-1]))
        shifts = [0.0]
        if np.any(ra0 + dra < 0.0):
            shifts.append(360.0)
        if np.any(ra0 + dra >= 360.0):
            shifts.append(-360.0)
        for s in shifts:
            lines.append(f"poly {len(verts)} {ra0 + s:.5f} {dec0:.5f} {jd:.5f} 1.0 500 {name}")
            lines += [f"  {x:+.5f} {y:+.5f}" for x, y in zip(dx, dy)]
            n_tiles += 1
    (HERE / "pointings.list").write_text("\n".join(lines) + "\n")
    return n_tiles, worst


def main(argv=None) -> int:
    with (HERE / "MPC.toml").open("rb") as fh:
        spec = tomllib.load(fh)["survey"]
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--h-lim", type=float, default=float(spec["h_lim"]))
    p.add_argument("--filter", default=spec["filter"])
    p.add_argument("--beta-max", type=float, default=8.0,
                   help="ecliptic latitude half-width [deg]; i_free < 4 deg orbits stay inside ~6.3")
    p.add_argument("--beta-sd", type=float, default=2.5,
                   help="latitude spread [deg] for the mean phase term")
    p.add_argument("--dlon", type=float, default=10.0)
    p.add_argument("--margin", type=float, default=1.0)
    p.add_argument("--width", type=float, default=0.05, help="tanh width [mag]")
    p.add_argument("--year-start", default="2026-01-01")
    args = p.parse_args(argv)

    n, worst = write_survey(args.h_lim, args.filter, args.beta_max, args.dlon, args.margin,
                            args.width, args.year_start, args.beta_sd)
    print(f"pointings.list: {n} tiles with their own .eff; max m0 error inside "
          f"37-62 au {worst:.3f} mag")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
