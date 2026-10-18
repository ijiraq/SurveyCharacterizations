"""Checks on the Eduardo et al. 2026 JWST Sample A characterization.

The layout tests are pure Python. The simulator tests need an ossssim that
loads multi-epoch survey roots (``read_survey_conf``); the line-of-sight
plant also needs ``survey_debias`` for its orbit helper.
"""
from __future__ import annotations

import csv
import math
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHAR = ROOT / "characterization"
NAME = "jwst_sample_a"
SURVEY = CHAR / NAME
# Copied verbatim from SurveySimulator-Data/Characterizations (single epoch).
SSDATA_SURVEYS = ("CFEPS", "OSSOS-MA", "All_Surveys", "All_r_Surveys", "Deep_Surveys")
EPOCH_JD = (2459969.32118, 2459973.96785, 2459979.90854)
FIELD_RA, FIELD_DEC = 209.3875, -10.865278
SIDE_DEG = math.sqrt(0.05)


def read_keyed(path: Path) -> dict[str, list[str]]:
    """``key= values`` lines of an .eff or survey.conf file (last one wins)."""
    out = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        out[key.strip()] = val.split()
    return out


def pointing_rows(path: Path) -> list[list[str]]:
    return [
        line.split() for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


class LayoutTest(unittest.TestCase):
    def test_survey_root_is_multi_epoch(self):
        self.assertFalse((SURVEY / "pointings.list").exists())
        epochs = sorted(p.name for p in SURVEY.iterdir()
                        if (p / "pointings.list").is_file())
        self.assertEqual(epochs, ["epoch1", "epoch2", "epoch3"])

    def test_survey_conf_requires_every_epoch(self):
        conf = read_keyed(SURVEY / "survey.conf")
        self.assertEqual(int(conf["detections_required"][0]), 3)

    def test_pointings(self):
        for i, jd in enumerate(EPOCH_JD, start=1):
            d = SURVEY / f"epoch{i}"
            rows = pointing_rows(d / "pointings.list")
            self.assertEqual(len(rows), 1)
            w, h, ra, dec, epoch, fill, observer, eff = rows[0]
            self.assertAlmostEqual(float(w), SIDE_DEG, places=6)
            self.assertAlmostEqual(float(h), SIDE_DEG, places=6)
            self.assertAlmostEqual(float(ra), FIELD_RA, places=6)
            self.assertAlmostEqual(float(dec), FIELD_DEC, places=6)
            self.assertAlmostEqual(float(epoch), jd, places=5)
            self.assertEqual(float(fill), 1.0)
            self.assertEqual(eff, f"{NAME}.eff")
            self.assertTrue((d / observer).is_file())
            self.assertTrue((d / eff).is_file())

    def test_observer_ephemeris_covers_epochs(self):
        for i, jd in enumerate(EPOCH_JD, start=1):
            text = (SURVEY / f"epoch{i}" / "JWST.csv").read_text()
            body = text.split("$$SOE")[1].split("$$EOE")[0]
            jds = [float(line.split(",")[0]) for line in body.splitlines() if line.strip()]
            self.assertLessEqual(min(jds), jd)
            self.assertGreaterEqual(max(jds), jd)

    def test_efficiency(self):
        effs = [(SURVEY / f"epoch{i}" / f"{NAME}.eff").read_text() for i in (1, 2, 3)]
        self.assertEqual(effs[0], effs[1])
        self.assertEqual(effs[0], effs[2])
        eff = read_keyed(SURVEY / "epoch1" / f"{NAME}.eff")
        self.assertEqual(eff["filter"], ["H"])
        self.assertEqual(eff["function"], ["single"])
        self.assertEqual([float(x) for x in eff["single_param"]], [0.96, 28.92, 0.61])
        self.assertEqual(float(eff["mag_lim"][0]), 30.5)
        self.assertEqual(float(eff["track_frac"][2]), 0.0)
        self.assertEqual([float(x) for x in eff["rate_cut"]], [0.03, 8.66, 0.0, 180.0])

    def test_mag_error_curve(self):
        e1, e2, e3, e4, _, _ = (float(x) for x in
                                read_keyed(SURVEY / "epoch1" / f"{NAME}.eff")["mag_error"])
        alpha = math.log10(e2 / e1) / (e3 - 21.0)

        def sigma(m):
            if m <= e3:
                return e1 * 10 ** (alpha * (m - 21.0))
            return max(e2 - (m - e3) * e4, 0.0)

        self.assertAlmostEqual(sigma(28.2), 0.05, places=3)
        self.assertAlmostEqual(sigma(29.5), 0.5, places=6)
        self.assertAlmostEqual(sigma(30.5), 0.5, places=6)

    def test_colour_maps_f150w2_to_r(self):
        with (ROOT / "colour.toml").open("rb") as fh:
            colours = tomllib.load(fh)
        for group in ("default", "cold"):
            self.assertAlmostEqual(colours[group]["r-g"] - colours[group]["H-g"], 1.0)

    def test_survey_config_and_detections(self):
        with (SURVEY / f"{NAME}.toml").open("rb") as fh:
            cfg = tomllib.load(fh)
        self.assertEqual(cfg["detection_defaults"], {"survey": NAME, "block": NAME})
        survey = cfg["survey"]
        self.assertEqual(survey["surveys"], [NAME])
        for key in ("detections_relpath", "results_name", "detections_full_name"):
            self.assertEqual(Path(survey[key]).parent, SURVEY.relative_to(ROOT), key)
        self.assertEqual(survey["detections_relpath"], f"characterization/{NAME}/{NAME}.csv")
        with (ROOT / survey["detections_relpath"]).open() as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 20)
        for col in cfg["detection_columns"].values():
            self.assertIn(col, rows[0])


class SurveySimulatorDataLayoutTest(unittest.TestCase):
    def test_single_epoch_survey_roots(self):
        for name in SSDATA_SURVEYS:
            with self.subTest(survey=name):
                d = CHAR / name
                rows = pointing_rows(d / "pointings.list")
                self.assertTrue(rows)
                self.assertFalse((d / "survey.conf").exists())
                effs = {row[-1] for row in rows if row[-1].endswith(".eff")}
                self.assertTrue(effs)
                for eff in effs:
                    self.assertTrue((d / eff).is_file(), eff)

    def test_filters_have_colours(self):
        with (ROOT / "colour.toml").open("rb") as fh:
            colours = tomllib.load(fh)
        for name in SSDATA_SURVEYS:
            for eff in (CHAR / name).glob("*.eff"):
                filt = read_keyed(eff).get("filter")
                if filt:
                    with self.subTest(eff=f"{name}/{eff.name}"):
                        self.assertIn(f"{filt[0]}-g", colours["default"])


def _ossssim_multi_epoch():
    try:
        import ossssim
    except ImportError:
        return None
    return ossssim if hasattr(ossssim, "read_survey_conf") else None


@unittest.skipIf(_ossssim_multi_epoch() is None,
                 "needs ossssim with multi-epoch survey roots")
class SimulatorTest(unittest.TestCase):
    def test_loads_as_three_epoch_survey(self):
        from ossssim import SurveyCharacterization

        s = SurveyCharacterization.from_directory(str(SURVEY))
        self.assertEqual(s.n_epochs, 3)
        self.assertEqual(s.detections_required, 3)
        self.assertEqual(len(s), 3)
        self.assertEqual({p.key for p in s}, {f"{NAME}/{NAME}"})

    def test_surveysimulator_data_loads_as_single_epoch(self):
        from ossssim import SurveyCharacterization

        for name in SSDATA_SURVEYS:
            with self.subTest(survey=name):
                s = SurveyCharacterization.from_directory(str(CHAR / name))
                self.assertEqual(s.n_epochs, 1)
                self.assertEqual(s.detections_required, 1)
                self.assertGreater(len(s), 0)
                self.assertEqual({p.survey for p in s}, {name})

    def test_line_of_sight_plant(self):
        try:
            from survey_debias.grid_bias import los_circular_elements
        except ImportError:
            self.skipTest("needs survey_debias")
        from astropy import units as u
        from ossssim import OSSSSim
        from ossssim.color import PhotSpec

        a, e, inc, node, peri, M = los_circular_elements(
            FIELD_RA, FIELD_DEC, 44.0, SURVEY / "epoch1" / "JWST.csv", EPOCH_JD[0]
        )
        sim = OSSSSim(str(SURVEY), seed=123)

        def flag(H):
            row = dict(a=a * u.au, e=e, inc=inc * u.deg, node=node * u.deg,
                       peri=peri * u.deg, M=M * u.deg, H=H * u.mag,
                       epoch=EPOCH_JD[0] * u.day, comp="default")
            out = sim.simulate(row, colors=PhotSpec(), model_band="r")
            return int(out["flag"]), str(out["Survey"]).strip()

        bright_flag, key = flag(8.0)
        self.assertGreaterEqual(bright_flag, 4)
        self.assertEqual(key, f"{NAME}/{NAME}")
        self.assertEqual(flag(16.0)[0], 0)


if __name__ == "__main__":
    unittest.main()
