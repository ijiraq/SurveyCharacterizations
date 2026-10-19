"""Write epoch{1,2,3}/pointings.list for jwst_sample_a from the NIRCam footprints.

footprints.csv holds one row per GO 1568 F150W2 SW calibrated exposure
(MAST service Mast.Jwst.Filtered.Nircam, productLevel 2b, *_cal.fits): the
detector's s_region and the exposure start (MJD). Observations 001, 002 and
003 are epochs 1, 2 and 3; visits are the 20 mosaic tiles; exposures 00001
and 00002 are the two dithers.

Usable sky, per tile and SW detector, is the part that detector imaged in both
dithers, less a 15-pixel edge band (Eduardo et al. 2026, Sections 4.3.3 and
5.5). Those 160 detector polygons per epoch overlap at tile edges, which
survey_debias's FoV-aimed sampling cannot handle, so each epoch is written as
a single polygon: the outline of the three NIRCam strips, joined across the
two 40" module gaps by 2" bridges. The fill factor is the fraction of that
outline on usable detector area.

    python make_pointings.py
"""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EFF = "jwst_sample_a.eff"
OBSERVER = "JWST.csv"
RA0, DEC0 = 209.3875, -10.865278          # Eduardo et al. field centre (ICRS)
EDGE_DEG = 15 * 0.031 / 3600.0            # 15 SW pixels
BRIDGE_DEG = 2.0 / 3600.0
EXPOSURE_DAYS = 639.0 / 86400.0           # one 3-integration exposure with resets
ROW_DEG = 0.3 / 60.0                      # tolerance when picking a strip's outer detector row
GRID_DEG = 1.0 / 3600.0                   # coverage sampling step
NAME = re.compile(r"jw01568(\d{3})(\d{3})_\d{5}_(\d{5})_(nrc[ab]\d)_cal\.fits")


def to_plane(ra, dec):
    """Gnomonic projection about (RA0, DEC0) [deg]."""
    a, d = np.radians(ra), np.radians(dec)
    a0, d0 = math.radians(RA0), math.radians(DEC0)
    c = math.sin(d0) * np.sin(d) + math.cos(d0) * np.cos(d) * np.cos(a - a0)
    x = np.cos(d) * np.sin(a - a0) / c
    y = (math.cos(d0) * np.sin(d) - math.sin(d0) * np.cos(d) * np.cos(a - a0)) / c
    return np.degrees(x), np.degrees(y)


def from_plane(x, y):
    x, y = np.radians(x), np.radians(y)
    a0, d0 = math.radians(RA0), math.radians(DEC0)
    rho = np.maximum(np.hypot(x, y), 1e-15)
    c = np.arctan(rho)
    d = np.arcsin(np.cos(c) * math.sin(d0) + y * np.sin(c) * math.cos(d0) / rho)
    a = a0 + np.arctan2(x * np.sin(c), rho * math.cos(d0) * np.cos(c) - y * math.sin(d0) * np.sin(c))
    return np.degrees(a) % 360.0, np.degrees(d)


def ccw(poly):
    x, y = poly[:, 0], poly[:, 1]
    signed = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    return poly if signed > 0 else poly[::-1]


def clip(poly, p, q, offset=0.0):
    """Part of convex poly left of the directed line p->q moved left by offset."""
    e = q - p
    n = np.array([-e[1], e[0]]) / np.hypot(*e)
    out = []
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        sa, sb = (a - p) @ n - offset, (b - p) @ n - offset
        if sa >= 0:
            out.append(a)
        if (sa >= 0) != (sb >= 0):
            out.append(a + sa / (sa - sb) * (b - a))
    return np.array(out)


def usable(dither1, dither2):
    """Sky imaged by one detector in both dithers, less the edge band."""
    both = ccw(dither1)
    for p, q in zip(ccw(dither2), np.roll(ccw(dither2), -1, axis=0)):
        both = clip(both, p, q)
    both = ccw(both)
    shrunk = both
    for p, q in zip(both, np.roll(both, -1, axis=0)):
        shrunk = clip(shrunk, p, q, EDGE_DEG)
    return ccw(shrunk)


def inside_convex(poly, pts):
    m = np.ones(len(pts), bool)
    for p, q in zip(poly, np.roll(poly, -1, axis=0)):
        e = q - p
        m &= (pts[:, 0] - p[0]) * e[1] - (pts[:, 1] - p[1]) * e[0] <= 0
    return m


def read_footprints(path):
    shots = defaultdict(dict)
    starts = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            m = NAME.fullmatch(row["filename"])
            if not m:
                continue
            obs, visit, expo, det = int(m[1]), int(m[2]), int(m[3]), m[4]
            v = np.array([float(t) for t in row["s_region"].split()[2:]]).reshape(-1, 2)
            shots[(obs, visit, det)][expo] = np.column_stack(to_plane(v[:, 0], v[:, 1]))
            starts[obs].append(float(row["expstart"]))
    return shots, starts


def strips(polys):
    """Three strips (bottom to top) as (axis, normal, along range, across range)."""
    allv = np.vstack(polys)
    w, v = np.linalg.eigh(np.cov(allv.T))
    u = v[:, np.argmax(w)]
    if u[0] < 0:
        u = -u
    nrm = np.array([-u[1], u[0]])
    cen = np.array([p.mean(0) @ nrm for p in polys])
    order = np.sort(cen)
    cuts = sorted(0.5 * (order[g] + order[g + 1]) for g in np.argsort(np.diff(order))[-2:])
    out = []
    for lo, hi in ((-np.inf, cuts[0]), (cuts[0], cuts[1]), (cuts[1], np.inf)):
        sel = [p for p, c in zip(polys, cen) if lo < c <= hi]
        along = np.concatenate([p @ u for p in sel])
        mins = np.array([np.min(p @ nrm) for p in sel])
        maxs = np.array([np.max(p @ nrm) for p in sel])
        out.append((along.min(), along.max(),
                    np.median(mins[mins < mins.min() + ROW_DEG]),
                    np.median(maxs[maxs > maxs.max() - ROW_DEG])))
    return u, nrm, sorted(out, key=lambda s: s[2])


def outline(u, nrm, s):
    """S-shaped simple polygon: bottom, middle and top strips, bridged at alternate ends."""
    (a00, a01, p00, p01), (a10, a11, p10, p11), (a20, a21, p20, p21) = s
    right = min(a01, a11)
    left = max(a10, a20)
    b = BRIDGE_DEG
    ap = [(a00, p00), (a01, p00), (a01, p01), (right, p01), (right, p10), (a11, p10),
          (a11, p11), (left + b, p11), (left + b, p20), (a21, p20), (a21, p21), (a20, p21),
          (a20, p20), (left, p20), (left, p11), (a10, p11), (a10, p10), (right - b, p10),
          (right - b, p01), (a00, p01)]
    pts = [a * u + p * nrm for a, p in ap]
    keep = [pts[0]]
    for q in pts[1:]:
        if np.hypot(*(q - keep[-1])) > 1e-9:
            keep.append(q)
    if np.hypot(*(keep[0] - keep[-1])) <= 1e-9:
        keep.pop()
    return np.array(keep)


def coverage(u, nrm, s, polys):
    """Fraction of the strip outline on usable detector area (grid sampling)."""
    covered = total = 0
    for a0, a1, p0, p1 in s:
        aa, pp = np.meshgrid(np.arange(a0, a1, GRID_DEG), np.arange(p0, p1, GRID_DEG))
        pts = aa.ravel()[:, None] * u + pp.ravel()[:, None] * nrm
        hit = np.zeros(len(pts), bool)
        for poly in polys:
            lo, hi = poly.min(0), poly.max(0)
            box = np.all((pts >= lo) & (pts <= hi), axis=1) & ~hit
            if box.any():
                hit[box] = inside_convex(poly, pts[box])
        covered += hit.sum()
        total += len(pts)
    return covered / total


def main():
    shots, starts = read_footprints(HERE / "footprints.csv")
    for obs in (1, 2, 3):
        polys = [usable(d[1], d[2]) for (o, _v, _d), d in sorted(shots.items()) if o == obs]
        u, nrm, s = strips(polys)
        poly = outline(u, nrm, s)
        fill = coverage(u, nrm, s, polys)
        x, y = poly[:, 0], poly[:, 1]
        area = 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)) * 3600.0
        jd = 0.5 * (min(starts[obs]) + max(starts[obs]) + EXPOSURE_DAYS) + 2400000.5
        ra_c, dec_c = from_plane(x.mean(), y.mean())
        ra, dec = from_plane(x, y)
        dx = (ra - ra_c) * np.cos(np.radians(dec))
        dy = dec - dec_c
        pa = math.degrees(math.atan2(u[1], u[0])) % 180.0
        widths = ", ".join(f"{(p1 - p0) * 60:.2f}" for _a0, _a1, p0, p1 in s)
        lines = [
            f"# JWST GO 1568 Sample A, epoch {obs} (observation {obs:03d}); written by ../make_pointings.py",
            f"# from ../footprints.csv. Outline of the three NIRCam SW strips ({widths}' wide, PA "
            f"{pa:.1f} deg in RA cos(dec), Dec),",
            f"# {area:.1f} arcmin^2 with the two module gaps excluded. Fill {fill:.2f}: {fill:.3f} of the "
            "outline was imaged by",
            "# one SW detector in both dithers, more than 15 pixels from its edge (README.md: a three-epoch",
            "# simulation of the detector polygons gives 0.913). JD: middle of the 20-tile epoch.",
            "# poly n ra dec jd fill observer efficiency; then n lines dx dy [deg], dx = dRA cos(dec)",
            f"poly {len(poly)} {ra_c:.7f} {dec_c:.7f} {jd:.5f} {fill:.2f} {OBSERVER} {EFF}",
        ]
        lines += [f"  {a:+.7f} {b:+.7f}" for a, b in zip(dx, dy)]
        out = HERE / f"epoch{obs}" / "pointings.list"
        out.write_text("\n".join(lines) + "\n")
        print(f"{out.relative_to(HERE)}: {len(poly)} vertices, {area:.1f} arcmin^2, "
              f"fill {fill:.3f}, JD {jd:.5f}")


if __name__ == "__main__":
    main()
