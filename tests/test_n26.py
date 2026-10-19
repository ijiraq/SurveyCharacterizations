"""Checks on the Napier et al. 2026 N26 HST heliostack characterization."""
from __future__ import annotations

import csv
import tomllib
import unittest
from pathlib import Path

from test_characterization import pointing_rows, read_keyed

ROOT = Path(__file__).resolve().parents[1]
NAME = "N26"
SURVEY = ROOT / "characterization" / NAME
EPOCH_JD = 2452672.8585


class N26LayoutTest(unittest.TestCase):
    def test_single_epoch_pointing(self):
        self.assertFalse((SURVEY / "survey.conf").exists())
        rows = pointing_rows(SURVEY / "pointings.list")
        self.assertEqual(len(rows), 1)
        shape, w, h, ra, dec, jd, fill, observer, eff = rows[0]
        self.assertEqual(shape, "rect")
        self.assertEqual((ra, dec), ("14:07:53.33", "-11:21:38"))
        self.assertAlmostEqual(float(jd), EPOCH_JD, places=4)
        self.assertEqual(float(fill), 1.0)
        self.assertEqual(eff, f"{NAME}.eff")
        self.assertTrue((SURVEY / observer).is_file())
        self.assertEqual({p.name for p in SURVEY.glob("*.eff")}, {eff})

    def test_observer_ephemeris_covers_epoch(self):
        body = (SURVEY / "HST.csv").read_text().split("$$SOE")[1].split("$$EOE")[0]
        jds = [float(line.split(",")[0]) for line in body.splitlines() if line.strip()]
        self.assertLessEqual(min(jds), EPOCH_JD)
        self.assertGreaterEqual(max(jds), EPOCH_JD)

    def test_efficiency_is_native_stmag(self):
        eff = read_keyed(SURVEY / f"{NAME}.eff")
        self.assertEqual(eff["filter"], ["t"])
        self.assertEqual(eff["function"], ["single"])
        self.assertEqual([float(x) for x in eff["single_param"]], [1.0, 29.21, 0.22])
        self.assertEqual([float(x) for x in eff["track_frac"]], [1.0, 29.21, -0.22])
        self.assertEqual(float(eff["mag_lim"][0]), 29.8)
        self.assertEqual([float(x) for x in eff["rate_cut"]], [0.05, 6.4, 0.0, 180.0])

    def test_colour_maps_f606w_stmag_to_r(self):
        with (ROOT / "colour.toml").open("rb") as fh:
            colours = tomllib.load(fh)
        for group in ("default", "cold"):
            self.assertAlmostEqual(colours[group]["t-g"] - colours[group]["r-g"], 0.3)

    def test_survey_config_and_detections(self):
        with (SURVEY / f"{NAME}.toml").open("rb") as fh:
            cfg = tomllib.load(fh)
        self.assertEqual(cfg["detection_defaults"], {"survey": NAME, "block": NAME})
        survey = cfg["survey"]
        self.assertEqual(survey["surveys"], [NAME])
        for key in ("detections_relpath", "results_name", "detections_full_name"):
            self.assertEqual(Path(survey[key]).parent, SURVEY.relative_to(ROOT), key)
        with (ROOT / survey["detections_relpath"]).open() as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual([r["name"] for r in rows],
                         ["2003BF91", "2003BH91", "2003ABCD", "2003WXYZ"])
        for col in cfg["detection_columns"].values():
            self.assertIn(col, rows[0])


if __name__ == "__main__":
    unittest.main()
