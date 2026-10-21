#!/usr/bin/env python
"""H_r of the MPC complete-survey objects from their reported photometry.

For each object (default: MPC H, h_mpc in MPC.csv, below --h-max and cold,
comp = cold, or a main-belt classical with a >= 42.4 au and i < 5 deg: the
sample and the objects near its limit under either cold definition) this

1. retrieves its observations from the MPC: the show_object page
   (https://minorplanetcenter.net/db_search/show_object?object_id=ID) creates
   the observation file (tmp2/ID.txt, 80-column format) and links to it; the
   file is cached in obs/<desig>.txt;
2. computes r, Delta and the phase angle alpha of every observation from the
   orbit in MPC.csv (two-body heliocentric propagation, light time; Earth from
   astropy) and checks the predicted positions against the observed ones;
3. reduces each magnitude in its reported band to
       H_X = m - 5 log10(r Delta) - beta alpha        (beta [mag/deg], linear)
   per photometric system: band@PS1 (F51, F52), band@LSST (X05), band@other;
4. moves H_X to H_r = H_X - (X - r), and takes the median of each night (per
   band, system and observatory), then the median over the nights.

The offsets X - r (``--offsets``):

  empirical (default)  measured from these objects: r@LSST is the reference;
                       a system band's offset is the median, over objects
                       observed in both, of median H_X - median H_r@LSST;
                       bands with no such object are tied to w@PS1 the same way.
  expected             cold colours from the root colour.toml rebased to r,
                       plus the bands it lacks (EXTRA_EXPECTED).

Only bands of the systems in ``--systems`` (default PS1 and LSST, calibrated
surveys) enter H_r; the others are listed. Writes MPC_photometry.csv,
MPC_band_offsets.txt and obs/<desig>_reduced.csv.

Usage (from this directory)::

    python mpc_photometry.py                       # cold objects with MPC H < 6
    python mpc_photometry.py --beta 0.10 --offsets expected
    python mpc_photometry.py --refresh             # re-download the observations
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
import tomllib
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OBS_DIR = HERE / "obs"
SHOW = "https://minorplanetcenter.net/db_search/show_object?object_id={}"
BASE = "https://minorplanetcenter.net/"
K_GAUSS = 0.01720209895                 # rad/day, heliocentric
OBLIQUITY = math.radians(23.4392911)    # J2000, the MPC's ecliptic
C_AU_DAY = 173.1446327
SYSTEMS = {"F51": "PS1", "F52": "PS1", "X05": "LSST"}
REFERENCE, BRIDGE = "r@LSST", "w@PS1"


def system_of(code: str) -> str:
    return SYSTEMS.get(code, "other")


def expected_offsets(colour_file: Path, group: str = "cold") -> dict:
    """X - r by band from colour.toml's group (rebased to r), plus EXTRA_EXPECTED."""
    with colour_file.open("rb") as fh:
        table = tomllib.load(fh)[group]
    vs_g = {k.split("-")[0]: float(v) for k, v in table.items()}
    off = {b: v - vs_g["r"] for b, v in vs_g.items()}
    gr, ri = off["g"], -off["i"]
    # ATLAS c ~ (g + r)/2, o ~ (r + i)/2 (Tonry et al. 2018); PS1 y ~ z - 0.05;
    # G: PS1 reports it with magnitudes on the r/w scale (see MPC_band_offsets.txt),
    # elsewhere Gaia G (Evans et al. 2018).
    off.setdefault("c", 0.5 * gr)
    off.setdefault("o", -0.5 * ri)
    off.setdefault("y", off.get("z", -ri) - 0.05)
    off.setdefault("G", -0.129 + 0.247 * gr - 0.027 * gr ** 2 - 0.049 * gr ** 3)
    return off


# ---------------------------------------------------------------------------
# Retrieval


def fetch(object_id: str, desig: str, refresh: bool = False, pause: float = 2.0) -> Path:
    """obs/<desig>.txt: the MPC observation file of ``object_id``."""
    OBS_DIR.mkdir(exist_ok=True)
    out = OBS_DIR / f"{desig}.txt"
    if out.is_file() and not refresh:
        return out
    headers = {"User-Agent": "SurveyCharacterizations/mpc_photometry"}
    url = SHOW.format(urllib.parse.quote(object_id))
    html = urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                  timeout=120).read().decode("utf-8", "replace")
    links = re.findall(r'(?:https?://[^"\s<>]*)?/?tmp2/[^"\s<>]+\.txt', html)
    if not links:
        raise RuntimeError(f"{object_id}: no tmp2 observation link on {url}")
    link = links[0] if links[0].startswith("http") else urllib.parse.urljoin(BASE, links[0].lstrip("/"))
    time.sleep(pause)
    text = urllib.request.urlopen(urllib.request.Request(link, headers=headers),
                                  timeout=120).read().decode("ascii", "replace")
    if not text.strip():
        raise RuntimeError(f"{object_id}: {link} is empty")
    out.write_text(f"# {link} retrieved {time.strftime('%Y-%m-%d')}\n" + text)
    time.sleep(pause)
    return out


# ---------------------------------------------------------------------------
# Observations and geometry


def read_obs(path: Path, kinds: str) -> list[dict]:
    """Observations with a magnitude from an 80-column MPC file (types in ``kinds``)."""
    rows = []
    for raw in path.read_text().splitlines():
        if raw.startswith("#") or len(raw) < 80 or raw[14] not in kinds:
            continue
        mag = raw[65:70].strip()
        if not mag:
            continue
        y, mo, d = raw[15:32].split()
        ra = raw[32:44].split()
        dec = raw[44:56].split()
        sign = -1.0 if dec[0].startswith("-") else 1.0
        rows.append({
            "date": f"{y}-{mo}-{d}", "y": int(y), "mo": int(mo), "day": float(d),
            "ra": 15.0 * (float(ra[0]) + float(ra[1]) / 60 + float(ra[2]) / 3600),
            "dec": sign * (abs(float(dec[0])) + float(dec[1]) / 60 + float(dec[2]) / 3600),
            "mag": float(mag), "band": raw[70].strip(), "code": raw[77:80], "kind": raw[14],
        })
    return rows


def jd_of(rows) -> np.ndarray:
    from astropy.time import Time

    iso = [f"{r['y']:04d}-{r['mo']:02d}-{int(r['day']):02d}" for r in rows]
    frac = np.array([r["day"] - int(r["day"]) for r in rows])
    return Time(iso, scale="utc").jd + frac


def helio_ecliptic(el: dict, jd: np.ndarray) -> np.ndarray:
    """Two-body heliocentric ecliptic J2000 position [au], shape (n, 3)."""
    a, e = el["a"], el["e"]
    inc, node, peri = (math.radians(el[k]) for k in ("i", "node", "peri"))
    M = math.radians(el["M"]) + K_GAUSS / a ** 1.5 * (jd - el["epoch_jd"])
    E = M.copy()
    for _ in range(50):
        E = E - (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
    x = a * (np.cos(E) - e)
    y = a * math.sqrt(1 - e * e) * np.sin(E)
    cw, sw, cn, sn, ci, si = (math.cos(peri), math.sin(peri), math.cos(node), math.sin(node),
                              math.cos(inc), math.sin(inc))
    px = (cw * cn - sw * sn * ci) * x + (-sw * cn - cw * sn * ci) * y
    py = (cw * sn + sw * cn * ci) * x + (-sw * sn + cw * cn * ci) * y
    pz = (sw * si) * x + (cw * si) * y
    return np.stack([px, py, pz], axis=1)


def ecl_to_eq(v: np.ndarray) -> np.ndarray:
    c, s = math.cos(OBLIQUITY), math.sin(OBLIQUITY)
    return np.stack([v[:, 0], c * v[:, 1] - s * v[:, 2], s * v[:, 1] + c * v[:, 2]], axis=1)


def geometry(el: dict, jd_utc: np.ndarray):
    """r, Delta [au], phase angle [deg], predicted RA, Dec [deg] (geocentric, light time)."""
    from astropy import units as u
    from astropy.coordinates import get_body_barycentric
    from astropy.time import Time

    t = Time(jd_utc, format="jd", scale="utc")
    earth = (get_body_barycentric("earth", t) - get_body_barycentric("sun", t)).xyz.to(u.au).value.T
    jd_tt = t.tt.jd
    obj = ecl_to_eq(helio_ecliptic(el, jd_tt))
    for _ in range(2):
        delta = np.linalg.norm(obj - earth, axis=1)
        obj = ecl_to_eq(helio_ecliptic(el, jd_tt - delta / C_AU_DAY))
    geo = obj - earth
    r = np.linalg.norm(obj, axis=1)
    delta = np.linalg.norm(geo, axis=1)
    alpha = np.degrees(np.arccos(np.clip(np.sum(obj * geo, axis=1) / (r * delta), -1.0, 1.0)))
    ra = np.degrees(np.arctan2(geo[:, 1], geo[:, 0])) % 360.0
    dec = np.degrees(np.arcsin(geo[:, 2] / delta))
    return r, delta, alpha, ra, dec


def sep_arcmin(ra1, dec1, ra2, dec2) -> np.ndarray:
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0))) * 60.0


# ---------------------------------------------------------------------------
# Reduction


def mad_sd(x) -> float:
    x = np.asarray(x, dtype=float)
    return 1.4826 * float(np.median(np.abs(x - np.median(x)))) if len(x) > 1 else float("nan")


def reduce_object(row: dict, path: Path, beta: float, kinds: str):
    """H_X of every observation with a band, and the orbit check."""
    from astropy.time import Time

    obs = [o for o in read_obs(path, kinds) if o["band"]]
    if not obs:
        return None, []
    el = {k: float(row[k]) for k in ("a", "e", "i", "node", "peri", "M")}
    el["epoch_jd"] = Time(row["epoch"], scale="tt").jd
    jd = jd_of(obs)
    r, delta, alpha, ra, dec = geometry(el, jd)
    offset = sep_arcmin(ra, dec, np.array([o["ra"] for o in obs]), np.array([o["dec"] for o in obs]))
    out = []
    for k, o in enumerate(obs):
        out.append({**o, "key": f"{o['band']}@{system_of(o['code'])}", "jd": jd[k],
                    "night": int(math.floor(jd[k] - 0.5)), "r": r[k], "delta": delta[k],
                    "alpha": alpha[k], "pos_offset_arcmin": offset[k],
                    "H_X": o["mag"] - 5.0 * math.log10(r[k] * delta[k]) - beta * alpha[k]})
    return {"max_offset": float(offset.max()), "median_offset": float(np.median(offset))}, out


def key_medians(obs: list[dict]) -> dict:
    """Median H_X per band@system over nights (night medians first)."""
    nights = {}
    for o in obs:
        nights.setdefault((o["key"], o["night"], o["code"]), []).append(o["H_X"])
    per_key = {}
    for (key, _, _), v in nights.items():
        per_key.setdefault(key, []).append(float(np.median(v)))
    return {k: float(np.median(v)) for k, v in per_key.items()}


def empirical_offsets(medians: dict) -> dict:
    """X - r per band@system: against r@LSST, else through w@PS1. (offset, sd, n objects)."""
    keys = sorted({k for m in medians.values() for k in m})

    def tie(key, ref):
        d = [m[key] - m[ref] for m in medians.values() if key in m and ref in m]
        return (float(np.median(d)), mad_sd(d), len(d)) if d else None

    out = {REFERENCE: (0.0, 0.0, sum(REFERENCE in m for m in medians.values()))}
    bridge = tie(BRIDGE, REFERENCE)
    for key in keys:
        if key == REFERENCE:
            continue
        direct = tie(key, REFERENCE)
        if direct is not None:
            out[key] = direct
        elif bridge is not None and key != BRIDGE:
            via = tie(key, BRIDGE)
            if via is not None:
                out[key] = (via[0] + bridge[0], via[1], via[2])
    return out


def mpc_h(row: dict) -> float:
    return float(row.get("h_mpc") or row["h"])


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--objects", nargs="*", default=None,
                   help="MPC.csv desig values (default: comp = cold and h_mpc < --h-max)")
    p.add_argument("--h-max", type=float, default=6.0, help="MPC H limit of the default list")
    p.add_argument("--beta", type=float, default=0.14,
                   help="linear phase coefficient [mag/deg] (TNOs: ~0.1-0.2)")
    p.add_argument("--offsets", choices=("empirical", "expected"), default="empirical")
    p.add_argument("--systems", default="PS1,LSST",
                   help="photometric systems used for H_r (PS1, LSST, other)")
    p.add_argument("--kinds", default="C", help="80-column observation types kept (C = CCD)")
    p.add_argument("--colour-file", default=str(ROOT / "colour.toml"))
    p.add_argument("--refresh", action="store_true", help="re-download the observations")
    p.add_argument("--out", default=str(HERE / "MPC_photometry.csv"))
    args = p.parse_args(argv)
    use = set(args.systems.split(","))

    lines = [x for x in (HERE / "MPC.csv").read_text().splitlines() if not x.startswith("#")]
    rows = {r["desig"]: r for r in csv.DictReader(lines)}
    targets = ([rows[d] for d in args.objects] if args.objects else
               [r for r in rows.values() if mpc_h(r) < args.h_max and (
                   r["comp"] == "cold" or (r["cl"] == "cla" and r["p"] == "m"
                                           and float(r["a"]) >= 42.4 and float(r["i"]) < 5.0))])
    targets.sort(key=mpc_h)

    reduced, infos = {}, {}
    for row in targets:
        object_id = row["number"] or re.sub(r"^(\d{4})", r"\1 ", row["prov"] or row["desig"])
        try:
            path = fetch(object_id, row["desig"], args.refresh)
        except (RuntimeError, OSError) as exc:
            print(f"{row['name']}: skipped ({exc})")
            continue
        info, obs = reduce_object(row, path, args.beta, args.kinds)
        if obs:
            reduced[row["desig"]], infos[row["desig"]] = obs, info

    medians = {d: key_medians(o) for d, o in reduced.items()}
    if args.offsets == "empirical":
        offsets = empirical_offsets(medians)
    else:
        exp = expected_offsets(Path(args.colour_file))
        offsets = {k: (exp[k.split("@")[0]], float("nan"), 0)
                   for m in medians.values() for k in m if k.split("@")[0] in exp}
    off_lines = [f"# X - r ({args.offsets}; phase beta = {args.beta:g} mag/deg); "
                 f"reference {REFERENCE}" + (f", bridge {BRIDGE}" if args.offsets == "empirical" else ""),
                 f"{'band@system':14s} {'X-r':>6s} {'sd':>5s} {'n_obj':>5s}"]
    for k, (v, sd, n) in sorted(offsets.items()):
        off_lines.append(f"{k:14s} {v:+6.2f} {sd:5.2f} {n:5d}")
    (HERE / "MPC_band_offsets.txt").write_text("\n".join(off_lines) + "\n")
    print("\n".join(off_lines) + "\n")

    summary = []
    for row in targets:
        desig = row["desig"]
        if desig not in reduced:
            print(f"{row['name']}: no usable photometry\n")
            continue
        obs, info = reduced[desig], infos[desig]
        for o in obs:
            v = offsets.get(o["key"])
            o["X_minus_r"] = v[0] if v else float("nan")
            o["H_r"] = o["H_X"] - o["X_minus_r"]
            o["used"] = int(v is not None and o["key"].split("@")[1] in use)
        write_reduced(obs, OBS_DIR / f"{desig}_reduced.csv")
        nights = {}
        for o in obs:
            nights.setdefault((o["key"], o["night"], o["code"]), []).append(o)
        per_night = [{"key": k[0], "code": k[2], "H_r": float(np.median([x["H_r"] for x in g])),
                      "used": g[0]["used"]} for k, g in nights.items()]
        used = np.array([n["H_r"] for n in per_night if n["used"]])
        print(f"{row['name']}: MPC H {mpc_h(row):.2f}; {len(obs)} observations; alpha "
              f"{min(o['alpha'] for o in obs):.2f}-{max(o['alpha'] for o in obs):.2f} deg; "
              f"orbit vs observed position: median {info['median_offset']:.1f}', "
              f"max {info['max_offset']:.1f}'")
        print(f"  {'band@system':14s} {'X-r':>6s} {'nights':>6s} {'H_r':>6s} {'sd':>5s}  codes")
        for key in sorted({n["key"] for n in per_night}):
            sel = [n for n in per_night if n["key"] == key]
            v = [n["H_r"] for n in sel]
            off = offsets.get(key, (float("nan"),))[0]
            flag = "" if sel[0]["used"] else "  (not used)"
            print(f"  {key:14s} {off:+6.2f} {len(sel):6d} {np.median(v):6.2f} {mad_sd(v):5.2f}  "
                  f"{','.join(sorted({n['code'] for n in sel}))}{flag}")
        if len(used) == 0:
            print("  no nights in the systems used\n")
            continue
        h_r = float(np.median(used))
        err = 1.2533 * mad_sd(used) / math.sqrt(len(used)) if len(used) > 1 else float("nan")
        print(f"  H_r = {h_r:.2f} +/- {err:.2f} ({len(used)} nights); MPC H - H_r = "
              f"{mpc_h(row) - h_r:+.2f}\n")
        summary.append({"desig": desig, "name": row["name"], "mpc_h": f"{mpc_h(row):.2f}",
                        "H_r": f"{h_r:.3f}", "H_r_err": f"{err:.3f}",
                        "sd_nights": f"{mad_sd(used):.3f}", "n_nights": len(used),
                        "bands": " ".join(sorted({n['key'] for n in per_night if n['used']})),
                        "alpha_max": f"{max(o['alpha'] for o in obs):.2f}"})
    with Path(args.out).open("w", newline="") as fh:
        fh.write(f"# H_r from MPC photometry (mpc_photometry.py): beta = {args.beta:g} mag/deg, "
                 f"{args.offsets} offsets (MPC_band_offsets.txt), systems {args.systems}\n")
        w = csv.DictWriter(fh, fieldnames=list(summary[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    diffs = np.array([float(s["mpc_h"]) - float(s["H_r"]) for s in summary])
    print(f"wrote {args.out}: {len(summary)} objects; MPC H - H_r median {np.median(diffs):+.2f}, "
          f"sd {np.std(diffs):.2f}")
    return 0


def write_reduced(rows, path: Path) -> None:
    cols = ["date", "code", "band", "key", "mag", "r", "delta", "alpha", "pos_offset_arcmin",
            "H_X", "X_minus_r", "H_r", "used"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.4f}" if isinstance(r[k], float) else r[k]) for k in cols})


if __name__ == "__main__":
    sys.exit(main())
