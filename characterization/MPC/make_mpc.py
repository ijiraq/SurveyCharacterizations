#!/usr/bin/env python
"""Build the MPC "complete survey" sample from MPCORB.DAT.

The Minor Planet Center orbit file is treated as a survey that has found
every object brighter than an absolute magnitude limit. This script writes

  MPC.csv          every MPCORB orbit with a_min <= a <= a_max and MPC H <= h_max,
                   with the sample columns used by survey_debias_joint; h is
                   H_r (see below), h_mpc the MPC's own H
  MPC.duplicates   MPC designation and the release (OSSOS, CFEPS, ...) name
                   of objects in both, for survey_debias_joint duplicates
  MPC_completeness.png / .txt (with --completeness) discovery year vs H

Usage (from this directory)::

    python make_mpc.py                         # MPCORB.DAT -> MPC.csv, MPC.duplicates
    python make_mpc.py --completeness          # MPC.csv -> MPC_completeness.{png,txt}
    python make_mpc.py --calibrate             # MPC.csv -> MPC_calibration.txt

MPCORB.DAT is not committed (about 320 MB); download it from
https://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT. Its snapshot
date (the latest last-observation date in the file) is recorded in the
MPC.csv header.

H columns: ``h`` is H_r on the PS1/Rubin r scale. Where MPC_photometry.csv
(mpc_photometry.py: the MPC observations reduced to H_r) has the object,
h is that measured H_r (``h_source`` = photometry, ``h_err`` its error);
otherwise h = h_mpc - offset, the offset being the median MPC H - H_r of
the measured objects (``h_source`` = mpc_offset, ``h_err`` the robust sd of
MPC H - H_r). A measured H_r more than REJECT from h_mpc - offset is not
used (``h_source`` = mpc_offset_photometry_rejected).

Free inclination: ``ifree`` is I_free of Huang, Gladman & Volk (2022, ApJS
259, 54; survey_debias.ifree): their Table 1 value where the object is
listed (``ifree_source`` = huang2022_table, ``huang_reso`` their dynamical
class: 0 classical, -1 scattering, pq resonant), else their forcing-pole
grid on the MPCORB elements made barycentric (huang2022_grid).

Sample columns: ``cl`` is "cla" for a non-resonant classical, "res" for a
resonant and "sca" for a scattering object, "" otherwise; ``p`` is the
OSSOS classical-belt part (i: a < 39.4, m: 39.4-47.7, o: a > 47.7).
``comp`` is "cold" where cl = cla, p = m, ifree < ifree_max (4 deg) and
a >= cold_a_min: the release cut detections_io.ColdCut(("cla",), ("m",),
None, ifree_max=4), restricted to the a range of the OSSOS 1.0 cold
components (42.4 au and out) so the sample matches the model the
efficiency is computed from. An object is "cla" if its a lies in the
main-belt range and q >= q_min, unless MPC.vetting (vet_resonances.py) or,
for the objects in their Table 1, Huang et al.'s own classification says
otherwise. ``secure`` is 1 for multi-opposition orbits;
the a, e, i of a one-opposition orbit (MPC uncertainty 9 is common) can
be far off, so near-limit insecure objects are listed in README.mpc.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHAR = HERE.parent
NCODE = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
CENTURY = {"I": 1800, "J": 1900, "K": 2000}
INNER_EDGE, OUTER_EDGE = 39.4, 47.7      # OSSOS classical belt parts (3:2 and 2:1)
COLD_A_MIN = 42.4                        # inner edge of the OSSOS 1.0 cold components
MPC_MINUS_R = 0.08                       # MPC H - H_r without MPC_photometry.csv
REJECT = 0.5                             # |MPC H - H_r - offset| above which photometry is not used
RELEASE_SURVEYS = ("OSSOS_v11", "CFEPS_r", "CFEPS_g", "CFEPS_presurvey", "MA", "HiLat")
COLUMNS = ("name", "desig", "number", "prov", "h", "h_err", "h_source", "h_mpc", "g", "a", "e",
           "i", "ifree", "ifree_source", "huang_reso", "q", "node", "peri",
           "M", "epoch", "nobs", "nopp", "first_year", "last_obs", "disc_year", "year_source",
           "secure", "cl", "p", "comp")
PROV = re.compile(r"^(\d{4}) ?([A-Z]{2}\d*)$")


def unpack_epoch(packed: str) -> str:
    """K2669 -> 2026-06-09."""
    year = CENTURY[packed[0]] + int(packed[1:3])
    return f"{year}-{NCODE.index(packed[3]):02d}-{NCODE.index(packed[4]):02d}"


def unpack_provisional(packed: str) -> str:
    """Packed provisional designation (K13SA0A) -> 2013SA100; numbers are returned as is."""
    packed = packed.strip()
    if re.fullmatch(r"\d+", packed):
        return str(int(packed))
    if len(packed) == 5 and packed[0] in NCODE[10:] and packed[1:].isdigit():
        return str(NCODE.index(packed[0]) * 10000 + int(packed[1:]))
    if len(packed) != 7 or packed[0] not in CENTURY:
        return packed
    year = CENTURY[packed[0]] + int(packed[1:3])
    cycle = NCODE.index(packed[4]) * 10 + int(packed[5])
    return f"{year}{packed[3]}{packed[6]}{cycle if cycle else ''}"


def parse_line(line: str) -> dict | None:
    """One MPCORB orbit (MPC export format), or None for header, blank or H-less lines."""
    if len(line) < 194 or not line[8:13].strip():
        return None
    try:
        a = float(line[92:103])
        e = float(line[70:79])
        inc = float(line[59:68])
    except ValueError:
        return None
    readable = line[166:194].strip()
    m = re.match(r"^\((\d+)\)\s*(.*)$", readable)
    number, rest = (m.group(1), m.group(2)) if m else ("", readable)
    pm = PROV.match(rest)
    prov = f"{pm.group(1)}{pm.group(2)}" if pm else ""
    arc = line[127:136].strip()
    last = line[194:202].strip()
    if "-" in arc:
        first_year = int(arc.split("-")[0])
    else:
        first_year = int(last[:4]) if last[:4].isdigit() else 0
    if prov:
        disc_year, source = int(prov[:4]), "provisional"
    else:
        disc_year, source = first_year, "first_obs"
    return {
        "name": readable, "desig": number or prov or line[0:7].strip(), "number": number,
        "prov": prov, "h": float(line[8:13]), "h_mpc": float(line[8:13]), "g": float(line[14:19] or 0.15),
        "a": a, "e": e, "i": inc, "q": a * (1.0 - e),
        "node": float(line[48:57]), "peri": float(line[37:46]), "M": float(line[26:35]),
        "epoch": unpack_epoch(line[20:25].strip()), "nobs": int(line[117:122].strip() or 0),
        "nopp": int(line[123:126].strip() or 0), "first_year": first_year, "last_obs": last,
        "disc_year": disc_year, "year_source": source,
        "secure": int(int(line[123:126].strip() or 0) >= 2),
    }


def read_mpcorb(path: Path, a_min: float, a_max: float, h_max: float):
    rows, snapshot = [], ""
    with path.open() as fh:
        for line in fh:
            try:
                if not a_min <= float(line[92:103]) <= a_max:
                    continue
            except ValueError:
                continue
            row = parse_line(line)
            if row is None or not a_min <= row["a"] <= a_max or row["h"] > h_max:
                continue
            snapshot = max(snapshot, row["last_obs"])
            rows.append(row)
    return rows, snapshot


def read_vetting(path: Path) -> dict:
    """desig -> (cl, note) from MPC.vetting (whitespace columns: desig cl note...)."""
    out = {}
    if not path.is_file():
        return out
    for raw in path.read_text().splitlines():
        words = raw.split("#", 1)[0].split(None, 2)
        if len(words) >= 2:
            out[words[0]] = (words[1], words[2] if len(words) > 2 else "")
    return out


def add_ifree(rows) -> None:
    """ifree, ifree_source, huang_reso for the main-belt rows (Huang et al. 2022).

    Huang et al. Table 1 (barycentric, double average) where the object is
    listed, by number or provisional designation; otherwise their forcing-pole
    grid on the MPCORB elements made barycentric at the MPCORB epoch.
    """
    import numpy as np
    from astropy.time import Time
    from survey_debias.ifree import free_inclination_helio, huang_lookup

    todo = []
    for r in rows:
        r["ifree"], r["ifree_source"], r["huang_reso"] = float("nan"), "", ""
        if not INNER_EDGE <= r["a"] <= OUTER_EDGE:
            continue
        hit = huang_lookup(r["number"], r["prov"])
        if hit is not None:
            r["huang_reso"] = hit["RESO"]
        if hit is not None and np.isfinite(hit["Ifree"]):
            r["ifree"], r["ifree_source"] = float(hit["Ifree"]), "huang2022_table"
        else:
            todo.append(r)
    if todo:
        jd = Time([r["epoch"] for r in todo], scale="tt").tdb.jd
        cols = {k: np.array([float(r[k]) for r in todo]) for k in ("a", "e", "i", "node", "peri", "M")}
        ifree = free_inclination_helio(cols["a"], cols["e"], cols["i"], cols["node"],
                                       cols["peri"], cols["M"], jd)
        for r, v in zip(todo, np.atleast_1d(ifree)):
            r["ifree"], r["ifree_source"] = float(v), "huang2022_grid"


def classify(rows, args, vetting: dict) -> None:
    for r in rows:
        a = r["a"]
        r["p"] = "i" if a < INNER_EDGE else ("m" if a <= OUTER_EDGE else "o")
        r["cl"] = "cla" if INNER_EDGE <= a <= OUTER_EDGE and r["q"] >= args.q_min else ""
        if r["desig"] in vetting:
            r["cl"] = vetting[r["desig"]][0]
        if r["cl"] and r.get("huang_reso") not in ("", None):
            # Huang et al. classified every object of their Table 1 by hand (10 Myr).
            r["cl"] = {0: "cla", -1: "sca"}.get(r["huang_reso"], "res")
        inc_ok = args.i_max is None or r["i"] < args.i_max
        free_ok = args.ifree_max is None or r["ifree"] < args.ifree_max
        r["comp"] = "cold" if (r["cl"] == "cla" and r["p"] == "m" and inc_ok and free_ok
                               and args.cold_a_min <= a <= args.cold_a_max) else ""


def release_rows(char_root: Path):
    """(survey, release row) for every release detection, as header -> value dicts."""
    for survey in RELEASE_SURVEYS:
        path = char_root / survey / f"{survey}.detections"
        if not path.is_file():
            continue
        header = None
        for raw in path.read_text().splitlines():
            if raw.startswith("#"):
                header = raw[1:].split()
                continue
            words = raw.split()
            if header is not None and len(words) == len(header):
                yield survey, dict(zip(header, words))


def release_names(char_root: Path) -> list[tuple[str, str, str]]:
    """(survey, object, MPC designation unpacked) for every release detection with an MPC id."""
    return [(survey, row["object"], unpack_provisional(row["MPC"]))
            for survey, row in release_rows(char_root) if row.get("MPC")]


def calibrate(rows, char_root: Path, out: Path) -> list[str]:
    """MPC H minus release Hsur per release filter, all objects and the cold sample."""
    import numpy as np

    by_key = {}
    for r in rows:
        for key in (r["prov"], r["number"]):
            if key:
                by_key[key] = r
    diffs = {}
    for survey, row in release_rows(char_root):
        mpc = unpack_provisional(row.get("MPC", ""))
        if mpc in by_key:
            r = by_key[mpc]
            hs = float(row["Hsur"])
            diffs.setdefault(row["Filt"], []).append((r["h_mpc"] - hs, r["h"] - hs,
                                                      r["comp"] == "cold"))
    lines = ["# MPC H (h_mpc) and MPC.csv H_r (h) minus release Hsur (OSSOS/CFEPS/MA/HiLat",
             "# .detections), by release filter: median, mean, std, robust sd",
             f"{'filt':>4s} {'sample':>6s} {'n':>3s}   {'h_mpc - Hsur':>26s}   {'h - Hsur':>26s}"]
    for filt, v in sorted(diffs.items()):
        cold = np.array([x[2] for x in v])
        for label, m in (("all", np.ones(len(v), bool)), ("cold", cold)):
            if m.sum() < 3:
                continue
            cells = []
            for col in (0, 1):
                x = np.array([y[col] for y in v])[m]
                med = float(np.median(x))
                rsd = 1.4826 * float(np.median(np.abs(x - med)))
                cells.append(f"{med:+6.2f} {x.mean():+6.2f} {x.std():5.2f} {rsd:5.2f}")
            lines.append(f"{filt:>4s} {label:>6s} {m.sum():3d}   {cells[0]:>26s}   {cells[1]:>26s}")
    out.write_text("\n".join(lines) + "\n")
    return lines


def write_duplicates(rows, char_root: Path, path: Path) -> int:
    by_key = {}
    for r in rows:
        for key in (r["prov"], r["number"]):
            if key:
                by_key[key] = r["desig"]
    groups = {}
    for survey, obj, mpc in release_names(char_root):
        if mpc in by_key:
            groups.setdefault(by_key[mpc], []).append(f"{obj}  # {survey}")
    with path.open("w") as fh:
        fh.write("# MPC designation and the release names of the same object\n"
                 "# (from the MPC column of each release .detections; make_mpc.py)\n")
        for desig in sorted(groups):
            for entry in groups[desig]:
                fh.write(f"{desig:12s} {entry}\n")
    return sum(len(v) for v in groups.values())


def to_r(rows, photometry: Path, reject: float = REJECT) -> dict:
    """Set h to H_r: measured (MPC_photometry.csv) where usable, else MPC H - offset.

    The offset is the median of MPC H - H_r over the measured objects; a
    measured H_r more than ``reject`` from MPC H - offset is not used.
    """
    import numpy as np

    measured = {}
    if photometry.is_file():
        lines = [x for x in photometry.read_text().splitlines() if not x.startswith("#")]
        measured = {r["desig"]: (float(r["H_r"]), float(r["H_r_err"]))
                    for r in csv.DictReader(lines)}
    by_desig = {r["desig"]: r for r in rows}
    diffs = np.array([by_desig[d]["h_mpc"] - h for d, (h, _) in measured.items() if d in by_desig])
    offset = float(np.median(diffs)) if len(diffs) else MPC_MINUS_R
    kept = np.abs(diffs - offset) <= reject
    if len(diffs):
        offset = float(np.median(diffs[kept]))
        scatter = 1.4826 * float(np.median(np.abs(diffs[kept] - offset)))
    else:
        scatter = float("nan")
    rejected = []
    for r in rows:
        r["h"], r["h_err"], r["h_source"] = r["h_mpc"] - offset, scatter, "mpc_offset"
        if r["desig"] in measured:
            h, err = measured[r["desig"]]
            if abs(r["h_mpc"] - h - offset) <= reject:
                r["h"], r["h_err"], r["h_source"] = h, err, "photometry"
            else:
                r["h_source"] = "mpc_offset_photometry_rejected"
                rejected.append(r["desig"])
    return {"offset": offset, "scatter": scatter, "n": int(kept.sum()) if len(diffs) else 0,
            "rejected": rejected}


def write_csv(rows, path: Path, meta: dict) -> None:
    with path.open("w", newline="") as fh:
        for k, v in meta.items():
            fh.write(f"# {k}: {v}\n")
        w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["h"], r["desig"])):
            out = dict(r)
            for k in ("a", "q"):
                out[k] = f"{r[k]:.4f}"
            for k in ("e",):
                out[k] = f"{r[k]:.5f}"
            for k in ("i", "node", "peri", "M"):
                out[k] = f"{r[k]:.4f}"
            out["h"], out["h_err"] = f"{r['h']:.3f}", f"{r['h_err']:.3f}"
            out["ifree"] = f"{r['ifree']:.4f}" if r["ifree"] == r["ifree"] else ""
            w.writerow(out)


def read_csv(path: Path) -> tuple[list[dict], dict]:
    meta, lines = {}, []
    for raw in path.read_text().splitlines():
        if raw.startswith("#"):
            k, _, v = raw[1:].partition(":")
            meta[k.strip()] = v.strip()
        elif raw.strip():
            lines.append(raw)
    rows = list(csv.DictReader(lines))
    for r in rows:
        for k in ("h", "h_mpc", "a", "e", "i", "q"):
            r[k] = float(r[k])
        r["ifree"] = float(r["ifree"]) if r.get("ifree") not in ("", None) else float("nan")
        r["disc_year"] = int(r["disc_year"])
    return rows, meta


# ---------------------------------------------------------------------------
# Completeness from discovery history


def completeness(rows, snapshot_year: int, recent: int, h_lo: float, h_hi: float,
                 dh: float, out_stem: Path) -> list[str]:
    """Cumulative discoveries vs year per H bin; the table goes to <stem>.txt."""
    import numpy as np

    edges = np.arange(h_lo, h_hi + 0.5 * dh, dh)
    since = snapshot_year - recent + 1
    sets = {"all (a 30-150 au)": rows, "cold (comp = cold)": [r for r in rows if r["comp"] == "cold"]}
    lines = [f"# MPC discoveries per H_r bin; 'recent' = discovered in {since}-{snapshot_year}",
             "# disc_year is the provisional-designation year (first observation for named",
             "# numbered objects, whose provisional designation MPCORB does not give)."]
    limits = {}
    for label, rs in sets.items():
        h = np.array([r["h"] for r in rs])
        y = np.array([r["disc_year"] for r in rs])
        lines.append(f"\n## {label}")
        lines.append(f"{'H_lo':>5s} {'H_hi':>5s} {'N':>5s} {'recent':>6s} {'last':>5s}")
        h_lim = h_lo
        clean = True
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (h >= lo) & (h < hi)
            n_recent = int((y[m] >= since).sum())
            last = int(y[m].max()) if m.any() else 0
            lines.append(f"{lo:5.2f} {hi:5.2f} {int(m.sum()):5d} {n_recent:6d} {last:5d}")
            if clean and n_recent == 0:
                h_lim = hi
            elif n_recent:
                clean = False
        limits[label] = h_lim
        lines.append(f"faintest H with no recent discovery in it or any brighter bin: {h_lim:.2f}")
    out_stem.with_suffix(".txt").write_text("\n".join(lines) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
    years = np.arange(1990, snapshot_year + 1)
    for ax, (label, rs) in zip(axes, sets.items()):
        h = np.array([r["h"] for r in rs])
        y = np.array([r["disc_year"] for r in rs])
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (h >= lo) & (h < hi)
            if m.sum() < 1:
                continue
            cum = np.array([(y[m] <= yr).sum() for yr in years]) / m.sum()
            ax.step(years, cum, where="post", label=f"{lo:.2f}-{hi:.2f} ({m.sum()})")
        ax.axvline(since - 0.5, color="k", ls=":", lw=0.8)
        ax.set_title(f"{label}: limit {limits[label]:.2f}")
        ax.set_xlabel("discovery year")
        ax.set_ylabel("fraction of today's sample known")
        ax.legend(fontsize=6, ncol=2, title="H_r")
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=130)
    return lines


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--mpcorb", default=str(HERE / "MPCORB.DAT"))
    p.add_argument("--out", default=str(HERE / "MPC.csv"))
    p.add_argument("--a-min", type=float, default=30.0, help="orbits kept in MPC.csv")
    p.add_argument("--a-max", type=float, default=150.0)
    p.add_argument("--h-max", type=float, default=7.0)
    p.add_argument("--q-min", type=float, default=37.0,
                   help="perihelion below which a main-belt object is not called classical")
    p.add_argument("--i-max", type=float, default=None, help="cold cut: ecliptic i < i_max")
    p.add_argument("--ifree-max", type=float, default=4.0,
                   help="cold cut: free inclination (Huang et al. 2022) < ifree_max")
    p.add_argument("--cold-a-min", type=float, default=COLD_A_MIN,
                   help="the OSSOS 1.0 cold components start at 42.4 au")
    p.add_argument("--cold-a-max", type=float, default=OUTER_EDGE)
    p.add_argument("--vetting", default=str(HERE / "MPC.vetting"))
    p.add_argument("--photometry", default=str(HERE / "MPC_photometry.csv"),
                   help="measured H_r (mpc_photometry.py) used for h")
    p.add_argument("--completeness", action="store_true",
                   help="read MPC.csv and write MPC_completeness.{txt,png}")
    p.add_argument("--recent", type=int, default=10, help="years counted as recent")
    p.add_argument("--dh", type=float, default=0.25)
    p.add_argument("--calibrate", action="store_true",
                   help="read MPC.csv and write MPC_calibration.txt (MPC H - release Hsur)")
    args = p.parse_args(argv)

    out = Path(args.out)
    if args.calibrate:
        rows, _ = read_csv(out)
        print("\n".join(calibrate(rows, CHAR, out.with_name("MPC_calibration.txt"))))
        return 0
    if args.completeness:
        import math
        import tomllib

        rows, meta = read_csv(out)
        snap = int(meta.get("mpcorb_snapshot", "2026")[:4])
        with (HERE / "MPC.toml").open("rb") as fh:
            h_lim = float(tomllib.load(fh)["survey"]["h_lim"])
        h_lo = h_lim - math.ceil((h_lim - 3.0) / args.dh) * args.dh   # h_lim is a bin edge
        lines = completeness(rows, snap, args.recent, h_lo, 7.0, args.dh,
                             out.with_name("MPC_completeness"))
        print("\n".join(lines))
        return 0

    rows, snapshot = read_mpcorb(Path(args.mpcorb), args.a_min, args.a_max, args.h_max)
    vetting = read_vetting(Path(args.vetting))
    add_ifree(rows)
    classify(rows, args, vetting)
    conv = to_r(rows, Path(args.photometry))
    meta = {
        "source": "MPCORB.DAT (https://www.minorplanetcenter.org/iau/MPCORB/MPCORB.DAT)",
        "mpcorb_snapshot": f"{snapshot[:4]}-{snapshot[4:6]}-{snapshot[6:8]}",
        "selection": f"{args.a_min:g} <= a <= {args.a_max:g} au, H <= {args.h_max:g}",
        "cold": (f"cl = cla (a {INNER_EDGE:g}-{OUTER_EDGE:g}, q >= {args.q_min:g}; Huang et al. "
                 "2022 classification where listed, else MPC.vetting), "
                 f"p = m, {args.cold_a_min:g} <= a <= {args.cold_a_max:g}"
                 + (f", i < {args.i_max:g} deg" if args.i_max is not None else "")
                 + (f", i_free < {args.ifree_max:g} deg" if args.ifree_max is not None else "")),
        "ifree": ("free inclination (Huang, Gladman & Volk 2022, ApJS 259, 54): their Table 1 "
                  "where listed (ifree_source = huang2022_table), else their forcing-pole grid "
                  "on the MPCORB elements made barycentric (huang2022_grid); main belt only"),
        "h": ("H_r (r band, PS1/Rubin r): measured from the MPC photometry where "
              "h_source = photometry (MPC_photometry.csv, mpc_photometry.py), else "
              f"h_mpc - {conv['offset']:.3f} (median MPC H - H_r of {conv['n']} measured "
              f"objects, robust sd {conv['scatter']:.3f} = h_err); photometry not used "
              f"(|MPC H - H_r - offset| > {REJECT:g}): {' '.join(conv['rejected']) or 'none'}"),
        "h_mpc": "MPC H as given in MPCORB (the MPC's own band corrections, close to V)",
        "i": "ecliptic J2000",
        "vetted": f"{len(vetting)} objects from {Path(args.vetting).name}",
    }
    write_csv(rows, out, meta)
    n_dup = write_duplicates(rows, CHAR, out.with_name("MPC.duplicates"))
    n_cold = sum(1 for r in rows if r["comp"] == "cold")
    sys.stderr.write(f"{out.name}: {len(rows)} orbits ({n_cold} cold), snapshot "
                     f"{meta['mpcorb_snapshot']}; MPC.duplicates: {n_dup} release names; "
                     f"h = H_r ({sum(r['h_source'] == 'photometry' for r in rows)} measured, "
                     f"others MPC H - {conv['offset']:.3f}, rejected {conv['rejected']})\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
