"""Write the DEEP B1 characterization: DEEP.csv and epoch{1..29}/{pointings.list,DEEP.eff}.

Inputs, all in source/:
  DEEP_VI_table2.dat, DEEP_VI_ReadMe   Smotherman et al. 2024 (AJ 167, 136) Table 2,
                                       VizieR J/AJ/167/136: the 110 linked TNOs
  DEEP_VI_table1.txt                   Smotherman et al. Table 1: the 29 long stares
  long_stares.csv                      NOIRLab Astro Data Archive, raw 120 s VR
                                       exposures of each long stare (fetch_long_stares)

The efficiency is Bernardinelli et al. 2024 (AJ 167, 134; "DEEP III") Eq. (3)
fitted to the survey-wide curve of their Figure 2 (fit_figure2), the same at
every night until the per-pointing-group parameters are available.

    python make_deep.py                 write DEEP.csv and the epoch directories
    python make_deep.py --fit-figure    refit Eq. (3) to the Figure 2 vector curve
    python make_deep.py --fetch         re-query the archive for source/long_stares.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "source"
NAME = "DEEP"
EFF = f"{NAME}.eff"
# Geocentric, as for OSSOS and CFEPS (ObsPos knows 500 or a JPL vector file,
# not 807): CTIO's parallax is 0.2" at 40 au, a few 0.01 "/hr in rate.
OBSERVER = "500"

# Reference epoch for d_bary: 2020-10-18.0, the middle of the 2020 run and of
# the 2019-2021 baseline over which m_VR is averaged.
REF_JD = 2459140.5
GAUSS_K = 0.01720209895
BARY_MASS = 1.00134                       # Sun + planets, for barycentric elements

# Cold classical cut, as in characterization/MPC/MPC.toml [cut].
COLD_A = (42.4, 47.7)
COLD_Q_MIN = 37.0
COLD_I_MAX = 5.0

# DEEP III Eq. (3), survey-wide fit (Figure 2, red curve): c / [(1 + exp(k1 x))(1 + exp(k2 x))],
# x = m - m25. OSSSSim's double tanh (A/4)(1 - tanh(x/d1))(1 - tanh(x/d2)) is
# identical with A = c, d = 2/k.
C_PEAK, M25, K1, K2 = 0.787, 26.22, 2.14, 5.02
# DEEP III Figure 2 (right): 50% rate completeness at 134.4 and 467.8 px/day.
PIXEL_ARCSEC = 0.263
R50 = (134.4, 467.8)
DIRECTION_HALF_WIDTH = 45.0               # |phi| < 45 deg about decreasing ecliptic longitude
MAG_LIM = 27.0
# Fit to the Table 2 uncertainties: sigma(m) = 0.10 * 10**(0.4 (m - 26)), 0.03 dex rms.
MAG_ERROR = (1.0e-3, 0.64, 28.0, 0.0, 40.0, 0.0)

# DECam CCD bounds in the gnomonic plane about the boresight [deg]
# (DESTNOSIM destnosim/des/ccd.py; N30 is dead and absent). DEEP searches
# every CCD except S30 and N30 (DESTNOSIM deep.py ccdmask = [2, 61]).
CCD_BOUNDS = {
    "N1": (-1.0811, -0.782681, -0.157306, -0.00750506),
    "N2": (-0.771362, -0.472493, -0.157385, -0.00749848),
    "N3": (-0.461205, -0.161464, -0.157448, -0.00749265),
    "N4": (-0.150127, 0.149894, -0.15747, -0.00749085),
    "N5": (0.161033, 0.460796, -0.157638, -0.0074294),
    "N6": (0.472171, 0.771045, -0.157286, -0.00740563),
    "N7": (0.782398, 1.08083, -0.157141, -0.0074798),
    "N8": (-0.92615, -0.627492, -0.321782, -0.172004),
    "N9": (-0.616455, -0.317043, -0.322077, -0.172189),
    "N10": (-0.305679, -0.00571999, -0.322071, -0.17217),
    "N11": (0.00565427, 0.305554, -0.322243, -0.172254),
    "N12": (0.31684, 0.616183, -0.322099, -0.172063),
    "N13": (0.627264, 0.925858, -0.321792, -0.171887),
    "N14": (-0.926057, -0.62726, -0.485961, -0.336213),
    "N15": (-0.616498, -0.317089, -0.486444, -0.336606),
    "N16": (-0.30558, -0.00578257, -0.486753, -0.336864),
    "N17": (0.00532179, 0.305123, -0.486814, -0.33687),
    "N18": (0.316662, 0.616018, -0.486495, -0.336537),
    "N19": (0.62708, 0.92578, -0.485992, -0.336061),
    "N20": (-0.770814, -0.471826, -0.650617, -0.500679),
    "N21": (-0.460777, -0.161224, -0.650817, -0.501097),
    "N22": (-0.149847, 0.149886, -0.650816, -0.501308),
    "N23": (0.161001, 0.460566, -0.650946, -0.501263),
    "N24": (0.47163, 0.770632, -0.650495, -0.500592),
    "N25": (-0.615548, -0.316352, -0.814774, -0.665052),
    "N26": (-0.305399, -0.00591217, -0.814862, -0.665489),
    "N27": (0.00550714, 0.304979, -0.815022, -0.665418),
    "N28": (0.316126, 0.615276, -0.814707, -0.664908),
    "N29": (-0.46018, -0.16101, -0.97887, -0.829315),
    "N31": (0.160884, 0.460147, -0.978775, -0.829426),
    "S1": (-1.08096, -0.782554, 0.00715956, 0.15689),
    "S2": (-0.7713, -0.47242, 0.0074194, 0.157269),
    "S3": (-0.4611, -0.161377, 0.00723009, 0.157192),
    "S4": (-0.149836, 0.150222, 0.00737069, 0.157441),
    "S5": (0.161297, 0.461031, 0.0072399, 0.1572),
    "S6": (0.472537, 0.771441, 0.00728934, 0.157137),
    "S7": (0.782516, 1.08097, 0.00742809, 0.15709),
    "S8": (-0.92583, -0.627259, 0.171786, 0.32173),
    "S9": (-0.616329, -0.31694, 0.171889, 0.321823),
    "S10": (-0.305695, -0.00579187, 0.172216, 0.322179),
    "S11": (0.00556739, 0.305472, 0.172237, 0.322278),
    "S12": (0.316973, 0.61631, 0.172015, 0.322057),
    "S13": (0.627389, 0.925972, 0.171749, 0.321672),
    "S14": (-0.925847, -0.627123, 0.335898, 0.48578),
    "S15": (-0.616201, -0.316839, 0.336498, 0.486438),
    "S16": (-0.305558, -0.00574858, 0.336904, 0.486749),
    "S17": (0.00557115, 0.305423, 0.33675, 0.486491),
    "S18": (0.316635, 0.615931, 0.33649, 0.486573),
    "S19": (0.627207, 0.925969, 0.336118, 0.485923),
    "S20": (-0.770675, -0.471718, 0.500411, 0.65042),
    "S21": (-0.46072, -0.161101, 0.501198, 0.650786),
    "S22": (-0.149915, 0.14982, 0.501334, 0.650856),
    "S23": (0.160973, 0.460482, 0.501075, 0.650896),
    "S24": (0.47167, 0.770647, 0.50045, 0.650441),
    "S25": (-0.615564, -0.316325, 0.66501, 0.814674),
    "S26": (-0.30512, -0.0056517, 0.665531, 0.81505),
    "S27": (0.00560886, 0.305082, 0.665509, 0.815022),
    "S28": (0.316158, 0.615391, 0.665058, 0.814732),
    "S29": (-0.46021, -0.160988, 0.829248, 0.978699),
    "S30": (-0.150043, 0.149464, 0.829007, 0.978648),
    "S31": (0.160898, 0.460111, 0.82932, 0.978804),
}
MASKED_CCDS = ("S30",)

# Table 2 numeric columns in order (DEEP_VI_ReadMe). 13 of the 110 rows are one
# byte short of the ReadMe layout, so rows are split on whitespace instead.
TABLE2 = ("a", "e", "i", "aop", "lan", "Tp", "e_a", "e_e", "e_i", "e_lan", "e_aop", "e_Tp",
          "VRmag", "VRMAG", "e_VRmag")
ARCHIVE = "https://astroarchive.noirlab.edu/api/adv_search/find/?limit=2000&rectype=file"
FIGURE2_PDF = HERE / "Bernardinelli_2024_AJ_167_134.pdf"


def jd_from_iso(stamp: str) -> float:
    t = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
    return 2440587.5 + t.timestamp() / 86400.0


def read_table2():
    rows = []
    for line in (SOURCE / "DEEP_VI_table2.dat").read_text().splitlines():
        words = line.split()
        if not words:
            continue
        row = {k: float(v) for k, v in zip(TABLE2, words)}
        flag = next(k for k, w in enumerate(words) if w in ("True", "False"))
        row["Class"] = " ".join(words[len(TABLE2):flag])
        row["Secure"] = words[flag] == "True"
        row["DEEP"] = words[flag + 1]
        row["MPC"] = " ".join(words[flag + 2:])
        rows.append(row)
    return rows


def bary_distance(a, e, tp_julian_year, jd):
    """Two-body barycentric distance [au] at jd from (a, e, Tp)."""
    tp_jd = 2451545.0 + (tp_julian_year - 2000.0) * 365.25
    n = GAUSS_K * math.sqrt(BARY_MASS) / a ** 1.5
    mean = math.fmod(n * (jd - tp_jd), 2.0 * math.pi)
    big_e = mean
    for _ in range(50):
        big_e -= (big_e - e * math.sin(big_e) - mean) / (1.0 - e * math.cos(big_e))
    return a * (1.0 - e * math.cos(big_e))


def composition(row) -> str:
    cls = row["Class"].split()[0].lower()
    if cls != "classical":
        return cls
    q = row["a"] * (1.0 - row["e"])
    cold = (COLD_A[0] <= row["a"] <= COLD_A[1] and q >= COLD_Q_MIN
            and row["i"] < COLD_I_MAX)
    return "cold" if cold else "hot"


def write_detections():
    # The Table 2 H_VR column (VRMAG) is not used: it is out of register with the
    # other columns (README.deep), so H is left to survey_debias (m_vr, d_bary).
    rows = read_table2()
    fields = ["name", "mpc", "m_vr", "m_err", "a", "e", "i", "Omega", "omega",
              "q", "d_bary", "class", "secure", "comp"]
    with (HERE / f"{NAME}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(rows, key=lambda r: r["VRmag"]):
            w.writerow([
                r["DEEP"], r["MPC"].replace(" ", ""), f"{r['VRmag']:.3f}", f"{r['e_VRmag']:.3f}",
                f"{r['a']:.4f}", f"{r['e']:.4f}", f"{r['i']:.3f}",
                f"{r['lan']:.3f}", f"{r['aop']:.3f}", f"{r['a'] * (1 - r['e']):.3f}",
                f"{bary_distance(r['a'], r['e'], r['Tp'], REF_JD):.3f}", r["Class"],
                str(r["Secure"]).lower(), composition(r),
            ])
    comps = {}
    for r in rows:
        comps[composition(r)] = comps.get(composition(r), 0) + 1
    print(f"{NAME}.csv: {len(rows)} detections, " + ", ".join(f"{v} {k}" for k, v in sorted(comps.items())))


def fetch_long_stares():
    """Query the NOIRLab archive for the raw VR exposures of every Table 1 long stare."""
    table1 = [line.split("\t") for line in (SOURCE / "DEEP_VI_table1.txt").read_text().splitlines()
              if line.startswith("B1")]
    out = []
    for field, night, ra, dec, *_ in table1:
        query = {
            "outfields": ["proposal", "dateobs_min", "dateobs_max", "ra_center", "dec_center", "ifilter"],
            "search": [["caldat", night, night], ["instrument", "decam"], ["proc_type", "raw"],
                       ["obs_type", "object"], ["exposure", 119, 121],
                       ["ra_center", float(ra) - 1.0, float(ra) + 1.0],
                       ["dec_center", float(dec) - 1.0, float(dec) + 1.0]],
        }
        req = urllib.request.Request(ARCHIVE, data=json.dumps(query).encode(),
                                     headers={"Content-Type": "application/json"})
        res = [r for r in json.load(urllib.request.urlopen(req, timeout=600))[1:]
               if r["ifilter"].startswith("VR")]
        res.sort(key=lambda r: r["dateobs_min"])
        out.append([field, night, len(res), res[0]["dateobs_min"], res[-1]["dateobs_max"],
                    f"{res[0]['ra_center']:.6f}", f"{res[0]['dec_center']:.6f}",
                    ";".join(sorted({r["proposal"] for r in res}))])
        print(*out[-1][:5])
    out.sort(key=lambda r: r[3])
    with (SOURCE / "long_stares.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["field", "night", "n_exp", "first_utc", "last_utc", "ra", "dec", "proposal"])
        w.writerows(out)


def decam_outline():
    """Stepped outline of the CCD rows, bridged across the inter-row gaps [deg]."""
    rows = {}
    for x0, x1, y0, y1 in CCD_BOUNDS.values():
        key = round(0.5 * (y0 + y1), 2)
        lo, hi, ylo, yhi = rows.get(key, (x0, x1, y0, y1))
        rows[key] = (min(lo, x0), max(hi, x1), min(ylo, y0), max(yhi, y1))
    rows = [rows[k] for k in sorted(rows)]
    cuts = [rows[0][2]] + [0.5 * (rows[k][3] + rows[k + 1][2]) for k in range(len(rows) - 1)] \
        + [rows[-1][3]]
    # Rows of equal length differ by < 0.001 deg at their ends; give them one edge.
    x_lo = [r[0] for r in rows]
    x_hi = [r[1] for r in rows]
    for k in range(1, len(rows)):
        if abs(x_lo[k] - x_lo[k - 1]) < 1e-3:
            x_lo[k] = x_lo[k - 1] = min(x_lo[k], x_lo[k - 1])
        if abs(x_hi[k] - x_hi[k - 1]) < 1e-3:
            x_hi[k] = x_hi[k - 1] = max(x_hi[k], x_hi[k - 1])
    right, left = [], []
    for k in range(len(rows)):
        right += [(x_hi[k], cuts[k]), (x_hi[k], cuts[k + 1])]
        left += [(x_lo[k], cuts[k]), (x_lo[k], cuts[k + 1])]
    pts = np.array(right + left[::-1])
    out = []
    for k in range(len(pts)):
        a, b, c = pts[k - 1], pts[k], pts[(k + 1) % len(pts)]
        if abs((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) > 1e-12:
            out.append(b)
    return np.array(out)


def polygon_area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def searched_area():
    return sum((x1 - x0) * (y1 - y0) for name, (x0, x1, y0, y1) in CCD_BOUNDS.items()
               if name not in MASKED_CCDS)


def ecliptic_west_angle(ra, dec):
    """Direction of decreasing ecliptic longitude at (ra, dec), in ossssim's
    rate_cut convention: degrees from due West, positive towards North."""
    eps = math.radians(23.4392911)

    def to_ecl(a, d):
        a, d = math.radians(a), math.radians(d)
        x, y, z = math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d)
        yy = y * math.cos(eps) + z * math.sin(eps)
        zz = -y * math.sin(eps) + z * math.cos(eps)
        return math.atan2(yy, x), math.asin(zz)

    def to_eq(lam, beta):
        x, y, z = math.cos(beta) * math.cos(lam), math.cos(beta) * math.sin(lam), math.sin(beta)
        yy = y * math.cos(eps) - z * math.sin(eps)
        zz = y * math.sin(eps) + z * math.cos(eps)
        return math.degrees(math.atan2(yy, x)) % 360.0, math.degrees(math.asin(zz))

    lam, beta = to_ecl(ra, dec)
    ra2, dec2 = to_eq(lam - math.radians(0.01), beta)
    west = ((ra - ra2 + 180.0) % 360.0 - 180.0) * math.cos(math.radians(dec))
    return math.degrees(math.atan2(dec2 - dec, west))


def eff_text(field, night, jd, angle):
    r_lo, r_hi = (r * PIXEL_ARCSEC / 24.0 for r in R50)
    d1, d2 = 2.0 / K1, 2.0 / K2
    return f"""\
# DEEP B1 long stare {field} {night} (JD {jd:.5f}); written by ../make_deep.py.
# Single-night selection of the KBMOD digital-tracking search (Smotherman et
# al. 2024, DEEP VI); survey.conf requires 4 nights for a linked orbit.
#
# PROVISIONAL: every night carries the survey-wide efficiency of Bernardinelli
# et al. 2024 (DEEP III) Figure 2. The per-pointing-group m25, k1, k2, c and the
# rate-function sharpness are in DEEP's pointing-group file, which is not
# public. Replace these lines night by night when it is available.
#
# Magnitudes are DECam VR calibrated to Pan-STARRS rP1 with no colour term
# (DEEP VI Section 2.1), labelled r.
#
# rate_cut= r_min r_max angle half_width
#   DEEP III Figure 2: 50% rate completeness at r50,1 = {R50[0]} and r50,2 = {R50[1]}
#   px/day (0.263 "/px), applied as a hard cut. The search grid is |phi| < 45 deg
#   about decreasing ecliptic longitude (DEEP III Section 2), {angle:+.1f} deg from
#   West (positive North) at this field centre.
#
rate_cut= {r_lo:.3f} {r_hi:.3f} {angle:.1f} {DIRECTION_HALF_WIDTH:.1f}
#
# mag_error: sigma(m) = 0.10 * 10**(0.4 (m - 26)), a fit to the Table 2 m_VR
# uncertainties (0.03 dex rms). Those are errors of the multi-night mean, so
# phot_frac draws one measurement.
#
mag_error= {' '.join(f'{v:g}' for v in MAG_ERROR)}
#
phot_frac= 1.0 0.0 0.0
#
# Linking is applied by survey.conf (4 nights); the 94% linking efficiency of
# DEEP VI Section 4.3 is one factor per object, not per night, and is not here.
#
track_frac= 1.0 {M25:.2f} 0.0
#
filter= r
#
rates=  0.00 20.00
#
# DEEP III Eq. (3): c / [(1 + exp(k1 (m - m25))) (1 + exp(k2 (m - m25)))],
# c = {C_PEAK}, m25 = {M25}, k1 = {K1}, k2 = {K2}: the survey-wide curve of Figure 2
# (vector path, residual rms 0.0002; c/2 at m = 25.96 as the paper quotes).
# double_param A M0 d1 d2 with A = c, M0 = m25, d = 2/k.
#
function= double
double_param= {C_PEAK:.3f} {M25:.2f} {d1:.4f} {d2:.4f}
#
# Implants reach m_VR = 28; the efficiency is below 0.003 at 27.0.
#
mag_lim= {MAG_LIM:.2f}
"""


def write_epochs():
    with (SOURCE / "long_stares.csv").open() as fh:
        stares = sorted(csv.DictReader(fh), key=lambda r: r["first_utc"])
    poly = decam_outline()
    fill = searched_area() / polygon_area(poly)
    for old in HERE.glob("epoch*"):
        shutil.rmtree(old)
    for k, s in enumerate(stares, start=1):
        jd = 0.5 * (jd_from_iso(s["first_utc"]) + jd_from_iso(s["last_utc"]))
        ra, dec = float(s["ra"]), float(s["dec"])
        angle = ecliptic_west_angle(ra, dec)
        out = HERE / f"epoch{k}"
        out.mkdir()
        lines = [
            f"# DEEP B1 long stare {k} of {len(stares)}: {s['field']} {s['night']}, {s['n_exp']} x 120 s VR; "
            "written by ../make_deep.py",
            "# from ../source/long_stares.csv. Boresight of the first exposure; JD: middle of the stare.",
            f"# Outline of the DECam CCD rows ({polygon_area(poly):.3f} deg^2), fill {fill:.3f}: the 60 searched",
            "# CCDs (every CCD but S30 and N30) cover that fraction of it.",
            "# poly n ra dec jd fill observer efficiency; then n lines dx dy [deg], dx = dRA cos(dec)",
            f"poly {len(poly)} {ra:.6f} {dec:.6f} {jd:.5f} {fill:.3f} {OBSERVER} {EFF}",
        ]
        lines += [f"  {x:+.6f} {y:+.6f}" for x, y in poly]
        (out / "pointings.list").write_text("\n".join(lines) + "\n")
        (out / EFF).write_text(eff_text(s["field"], s["night"], jd, angle))
    print(f"epoch1..epoch{len(stares)}: DECam outline {len(poly)} vertices, "
          f"{polygon_area(poly):.3f} deg^2, fill {fill:.3f} ({searched_area():.3f} deg^2 searched)")


def fit_figure2():
    """Fit DEEP III Eq. (3) to the survey-wide (red) curve of Figure 2."""
    import pymupdf
    from scipy.optimize import least_squares

    page = pymupdf.open(FIGURE2_PDF)[4]
    curve = next(d for d in page.get_drawings()
                 if d.get("color") and np.allclose(d["color"], (0.74, 0.0, 0.1), atol=0.01)
                 and d["type"] == "s" and len(d["items"]) > 10)
    # Axes: m 22..27 at x = 90.459..293.371 pt; completeness 0..1 at y = 216.06..55.67 pt.
    pts = []
    for item in curve["items"]:
        ctrl = [np.array([p.x, p.y]) for p in item[1:]]
        if item[0] == "l":
            pts += ctrl
        else:
            for t in np.linspace(0.0, 1.0, 40):
                pts.append((1 - t) ** 3 * ctrl[0] + 3 * (1 - t) ** 2 * t * ctrl[1]
                           + 3 * (1 - t) * t ** 2 * ctrl[2] + t ** 3 * ctrl[3])
    pts = np.array(pts)
    m = 22.0 + 5.0 * (pts[:, 0] - 90.459) / (293.371 - 90.459)
    eta = (216.06 - pts[:, 1]) / (216.06 - 55.67)

    def resid(p):
        c, k1, k2 = p
        return c / ((1 + np.exp(k1 * (m - M25))) * (1 + np.exp(k2 * (m - M25)))) - eta

    fit = least_squares(resid, [0.8, 3.0, 10.0])
    c, k1, k2 = fit.x
    print(f"c = {c:.3f}, k1 = {k1:.2f}, k2 = {k2:.2f} (m25 = {M25}); "
          f"rms {np.sqrt(np.mean(fit.fun ** 2)):.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fetch", action="store_true", help="re-query the archive first")
    ap.add_argument("--fit-figure", action="store_true", help="refit Eq. (3) to Figure 2 and stop")
    args = ap.parse_args()
    if args.fit_figure:
        fit_figure2()
        return
    if args.fetch:
        fetch_long_stares()
    write_detections()
    write_epochs()


if __name__ == "__main__":
    main()
