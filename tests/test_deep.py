"""Checks on the DEEP B1 linked-sample characterization (Smotherman et al. 2024)."""
from __future__ import annotations

import csv
import importlib.util
import math
import tomllib
import unittest
from pathlib import Path

from test_characterization import poly_pointings, read_keyed

ROOT = Path(__file__).resolve().parents[1]
NAME = "DEEP"
SURVEY = ROOT / "characterization" / NAME
N_EPOCHS = 29
# DEEP III Eq. (3) fitted to Figure 2: c, m25, k1, k2.
C, M25, K1, K2 = 0.787, 26.22, 2.14, 5.02
COMPS = {"cold": 73, "hot": 23, "resonant": 6, "detached": 6, "scattering": 2}


def load_make_deep():
    spec = importlib.util.spec_from_file_location("make_deep", SURVEY / "make_deep.py")
    md = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(md)
    return md


def stares():
    with (SURVEY / "source" / "long_stares.csv").open() as fh:
        return sorted(csv.DictReader(fh), key=lambda r: r["first_utc"])


class DeepLayoutTest(unittest.TestCase):
    def test_multi_epoch_root(self):
        self.assertFalse((SURVEY / "pointings.list").exists())
        dirs = sorted(p.name for p in SURVEY.iterdir() if (p / "pointings.list").is_file())
        self.assertEqual(sorted(dirs), sorted(f"epoch{i}" for i in range(1, N_EPOCHS + 1)))
        self.assertEqual(read_keyed(SURVEY / "survey.conf")["detections_required"], ["4"])

    def test_pointings_follow_the_long_stares(self):
        md = load_make_deep()
        jds = []
        for k, stare in enumerate(stares(), start=1):
            with self.subTest(epoch=k):
                d = SURVEY / f"epoch{k}"
                polys = poly_pointings(d / "pointings.list")
                self.assertEqual(len(polys), 1)
                head, verts = polys[0]
                _, n, ra, dec, jd, fill, observer, eff = head
                self.assertEqual(int(n), len(verts))
                self.assertAlmostEqual(float(ra), float(stare["ra"]), places=5)
                self.assertAlmostEqual(float(dec), float(stare["dec"]), places=5)
                mid = 0.5 * (md.jd_from_iso(stare["first_utc"]) + md.jd_from_iso(stare["last_utc"]))
                self.assertAlmostEqual(float(jd), mid, places=4)
                self.assertEqual((fill, observer, eff), ("0.860", "500", f"{NAME}.eff"))
                self.assertEqual({p.name for p in d.glob("*.eff")}, {eff})
                jds.append(float(jd))
        self.assertEqual(jds, sorted(jds))

    def test_stares_match_deep_vi_table1(self):
        table1 = {}
        for line in (SURVEY / "source" / "DEEP_VI_table1.txt").read_text().splitlines():
            if line.startswith("B1"):
                field, night, ra, dec, n = line.split("\t")[:5]
                table1[(field, night)] = (float(ra), float(dec), int(n))
        self.assertEqual(len(table1), N_EPOCHS)
        for s in stares():
            with self.subTest(stare=(s["field"], s["night"])):
                ra, dec, n = table1[(s["field"], s["night"])]
                sep = math.hypot((float(s["ra"]) - ra) * math.cos(math.radians(dec)),
                                 float(s["dec"]) - dec) * 3600
                self.assertLess(sep, 10.0)
                self.assertLessEqual(abs(int(s["n_exp"]) - n), 1)

    def test_footprint_is_the_decam_outline(self):
        md = load_make_deep()
        poly = md.decam_outline()
        self.assertAlmostEqual(md.polygon_area(poly), 3.130, places=3)
        self.assertAlmostEqual(md.searched_area() / md.polygon_area(poly), 0.860, places=3)
        _, verts = poly_pointings(SURVEY / "epoch1" / "pointings.list")[0]
        self.assertEqual(len(verts), len(poly))
        for (x, y), (px, py) in zip(verts, poly):
            self.assertAlmostEqual(x, px, places=5)
            self.assertAlmostEqual(y, py, places=5)

    def test_efficiency_is_deep_iii_eq3(self):
        for k in range(1, N_EPOCHS + 1):
            eff = read_keyed(SURVEY / f"epoch{k}" / f"{NAME}.eff")
            with self.subTest(epoch=k):
                self.assertEqual(eff["filter"], ["r"])
                self.assertEqual(eff["function"], ["double"])
                a, m0, d1, d2 = (float(x) for x in eff["double_param"])
                self.assertEqual((a, m0), (C, M25))
                self.assertAlmostEqual(d1, 2 / K1, places=4)
                self.assertAlmostEqual(d2, 2 / K2, places=4)
                r_lo, r_hi, angle, width = (float(x) for x in eff["rate_cut"])
                self.assertAlmostEqual(r_lo, 134.4 * 0.263 / 24, places=3)
                self.assertAlmostEqual(r_hi, 467.8 * 0.263 / 24, places=3)
                self.assertAlmostEqual(angle, -23.3, delta=0.3)
                self.assertEqual(width, 45.0)
                self.assertEqual(float(eff["mag_lim"][0]), 27.0)
                self.assertEqual([float(x) for x in eff["phot_frac"]], [1.0, 0.0, 0.0])

    def test_double_tanh_reproduces_the_paper(self):
        def eq3(m):
            x = m - M25
            return C / ((1 + math.exp(K1 * x)) * (1 + math.exp(K2 * x)))

        def double(m):
            x = m - M25
            return C / 4 * (1 - math.tanh(x / (2 / K1))) * (1 - math.tanh(x / (2 / K2)))

        for m in (22.0, 25.0, 25.96, 26.22, 26.5, 27.0):
            self.assertAlmostEqual(eq3(m), double(m), places=12)
        self.assertAlmostEqual(eq3(M25), C / 4, places=12)
        self.assertAlmostEqual(eq3(25.96), C / 2, delta=0.002)   # DEEP III's converted m50
        self.assertLess(eq3(27.0), 0.003)

    def test_detections(self):
        with (SURVEY / f"{NAME}.toml").open("rb") as fh:
            cfg = tomllib.load(fh)
        self.assertEqual(cfg["detection_defaults"], {"survey": NAME, "block": NAME})
        survey = cfg["survey"]
        self.assertEqual(survey["surveys"], [NAME])
        self.assertEqual(survey["bias_method"], "aq_grid")
        for key in ("detections_relpath", "results_name", "detections_full_name"):
            self.assertEqual(Path(survey[key]).parent, SURVEY.relative_to(ROOT), key)
        with (ROOT / survey["detections_relpath"]).open() as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 110)
        self.assertEqual(len({r["name"] for r in rows}), 110)
        for col in cfg["detection_columns"].values():
            self.assertIn(col, rows[0])
        self.assertNotIn("h_vr", rows[0])
        counts = {}
        for r in rows:
            counts[r["comp"]] = counts.get(r["comp"], 0) + 1
            if r["comp"] == "cold":
                self.assertLess(float(r["i"]), 5.0)
                self.assertTrue(42.4 <= float(r["a"]) <= 47.7)
                self.assertGreaterEqual(float(r["q"]), 37.0)
        self.assertEqual(counts, COMPS)

    def test_magnitudes_give_jpl_h(self):
        # m_vr - 5 log10(r (r - 1)) at d_bary against JPL SBDB H (README.deep):
        # this is why DEEP VI's H_VR column is not used.
        jpl = {"2019QC110": 7.72, "2021RW237": 6.20, "2020UT74": 9.00, "2002PV170": 6.20}
        with (SURVEY / f"{NAME}.csv").open() as fh:
            rows = {r["mpc"]: r for r in csv.DictReader(fh)}
        for mpc, h in jpl.items():
            r = rows[mpc]
            d = float(r["d_bary"])
            with self.subTest(mpc=mpc):
                self.assertAlmostEqual(float(r["m_vr"]) - 5 * math.log10(d * (d - 0.99)), h, delta=0.15)


def _ossssim_multi_epoch():
    try:
        import ossssim
    except ImportError:
        return None
    return ossssim if hasattr(ossssim, "read_survey_conf") else None


@unittest.skipIf(_ossssim_multi_epoch() is None, "needs ossssim with multi-epoch survey roots")
class DeepSimulatorTest(unittest.TestCase):
    def test_loads_as_29_epoch_survey(self):
        from ossssim import SurveyCharacterization

        s = SurveyCharacterization.from_directory(str(SURVEY))
        self.assertEqual(s.n_epochs, N_EPOCHS)
        self.assertEqual(s.detections_required, 4)
        self.assertEqual(len(s), N_EPOCHS)
        self.assertEqual({p.key for p in s}, {f"{NAME}/{NAME}"})


if __name__ == "__main__":
    unittest.main()
