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
# Single-epoch characterizations from the OSSOS v11 release, one survey each.
SINGLE_EPOCH_SURVEYS = ("CFEPS_g", "CFEPS_r", "CFEPS_presurvey", "HiLat", "MA", "OSSOS_v11")
EPOCH_JD = (2459969.32121, 2459973.96788, 2459979.90857)
FIELD_RA, FIELD_DEC = 209.3875, -10.865278
FILL = 0.91
# Eduardo et al. 2026 Eq. (8): single-epoch p_det = p0/2 erfc((m - m0)/w).
P0, M0, W = 0.96, 28.92, 0.61
# Table 1 sample A; JPB08 is cold (38 < d < 55 au, i < 5 deg), JPB09 (d < 30 au)
# is the centaur left out of the luminosity-function samples.
COLD = {"JPB01", "JPB04", "JPB05", "JPB06", "JPB07", "JPB08", "JPB12", "JPB13",
        "JPB15", "JPB16", "JPB20", "JPB21"}
HOT = {"JPB10", "JPB11", "JPB14", "JPB17", "JPB18", "JPB19", "JPB22"}
CENTAUR = {"JPB09"}


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


def poly_pointings(path: Path) -> list[tuple[list[str], list[tuple[float, float]]]]:
    """``poly`` lines of a pointings.list with their (dx, dy) vertex lines."""
    rows = pointing_rows(path)
    out = []
    k = 0
    while k < len(rows):
        head = rows[k]
        n = int(head[1])
        out.append((head, [(float(a), float(b)) for a, b in rows[k + 1:k + 1 + n]]))
        k += 1 + n
    return out


def lookup_table(path: Path) -> list[tuple[float, float]]:
    return [
        tuple(float(v) for v in line.partition("=")[2].split()[:2])
        for line in path.read_text().splitlines() if line.startswith("lookup_param=")
    ]


def detection_objects(path: Path) -> list[str]:
    """Object column of a CFEPS-format .detections file."""
    return [row[5] for row in pointing_rows(path) if len(row) > 20]


def cds_objects(path: Path) -> list[str]:
    """Object column (bytes 18-27) of the data rows of a .CDS file."""
    return [
        line[17:27].strip() for line in path.read_text().splitlines()
        if len(line) > 250 and line[:3].isalpha() and line[:3].islower() and line[3] == " "
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
            polys = poly_pointings(d / "pointings.list")
            # One polygon per epoch: survey_debias aims plants per tile and
            # needs tiles that do not overlap.
            self.assertEqual(len(polys), 1)
            head, vertices = polys[0]
            _, n, ra, dec, epoch, fill, observer, eff = head
            self.assertEqual(int(n), len(vertices))
            self.assertLess(len(vertices), 41)
            self.assertAlmostEqual(float(ra), FIELD_RA, delta=0.01)
            self.assertAlmostEqual(float(dec), FIELD_DEC, delta=0.01)
            self.assertAlmostEqual(float(epoch), jd, places=5)
            self.assertEqual(float(fill), FILL)
            self.assertEqual(eff, f"{NAME}.eff")
            self.assertTrue((d / observer).is_file())
            self.assertTrue((d / eff).is_file())
            x = [v[0] for v in vertices]
            y = [v[1] for v in vertices]
            area = 0.5 * abs(sum(x[k] * y[k - 1] - x[k - 1] * y[k] for k in range(len(x))))
            # Three strips, 2.19' + 4.16' + 2.19' by 20.9', module gaps excluded.
            self.assertAlmostEqual(area * 3600, 178.9, delta=0.5)

    def test_footprints_cover_every_exposure(self):
        with (SURVEY / "footprints.csv").open() as fh:
            rows = list(csv.DictReader(fh))
        # 3 epochs x 20 tiles x 2 dithers x 8 SW detectors.
        self.assertEqual(len(rows), 960)
        self.assertEqual({int(r["observtn"]) for r in rows}, {1, 2, 3})
        self.assertEqual({int(r["visit"]) for r in rows}, set(range(1, 21)))

    def test_make_pointings_reproduces_epoch1(self):
        import importlib.util

        spec = importlib.util.spec_from_file_location("make_pointings", SURVEY / "make_pointings.py")
        mp = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mp)
        shots, _ = mp.read_footprints(SURVEY / "footprints.csv")
        polys = [mp.usable(d[1], d[2]) for (o, _v, _d), d in sorted(shots.items()) if o == 1]
        self.assertEqual(len(polys), 160)
        u, nrm, strips = mp.strips(polys)
        poly = mp.outline(u, nrm, strips)
        self.assertAlmostEqual(mp.coverage(u, nrm, strips, polys), 0.905, delta=0.002)
        head, vertices = poly_pointings(SURVEY / "epoch1" / "pointings.list")[0]
        ra, dec = mp.from_plane(poly[:, 0], poly[:, 1])
        dec0 = float(head[3])
        ra0 = float(head[2])
        for (dx, dy), r, d in zip(vertices, ra, dec):
            self.assertAlmostEqual(dec0 + dy, d, places=6)
            self.assertAlmostEqual(ra0 + dx / math.cos(math.radians(dec0 + dy)), r, places=6)

    def test_field_centre_is_imaged(self):
        """Eduardo et al.'s field centre lies in the middle strip of every epoch."""
        for i in (1, 2, 3):
            head, vertices = poly_pointings(SURVEY / f"epoch{i}" / "pointings.list")[0]
            ra0, dec0 = float(head[2]), float(head[3])
            px = (FIELD_RA - ra0) * math.cos(math.radians(FIELD_DEC))
            py = FIELD_DEC - dec0
            inside = False
            for (x1, y1), (x2, y2) in zip(vertices, vertices[1:] + vertices[:1]):
                if (y1 > py) != (y2 > py) and px < x1 + (py - y1) * (x2 - x1) / (y2 - y1):
                    inside = not inside
            self.assertTrue(inside, f"epoch{i}")

    def test_observer_ephemeris_covers_epochs(self):
        for i, jd in enumerate(EPOCH_JD, start=1):
            text = (SURVEY / f"epoch{i}" / "JWST.csv").read_text()
            body = text.split("$$SOE")[1].split("$$EOE")[0]
            jds = [float(line.split(",")[0]) for line in body.splitlines() if line.strip()]
            self.assertLessEqual(min(jds), jd)
            self.assertGreaterEqual(max(jds), jd)

    def test_efficiency(self):
        """Same single-epoch efficiency at every epoch; only the KBMOD rate grid differs."""
        texts = []
        for i in (1, 2, 3):
            path = SURVEY / f"epoch{i}" / f"{NAME}.eff"
            texts.append([line for line in path.read_text().splitlines()
                          if "=" in line and not line.startswith("#")
                          and not line.startswith("rate_cut=")])
            eff = read_keyed(path)
            self.assertEqual(eff["filter"], ["H"])
            self.assertEqual(eff["function"], ["lookup"])
            self.assertEqual(float(eff["mag_lim"][0]), 30.5)
            self.assertEqual(float(eff["track_frac"][2]), 0.0)
            r_max = 2.8 if i == 1 else 1.4
            self.assertEqual([float(x) for x in eff["rate_cut"]], [0.03, r_max, 0.0, 180.0])
        self.assertEqual(texts[0], texts[1])
        self.assertEqual(texts[0], texts[2])

    def test_lookup_is_paper_erfc(self):
        table = lookup_table(SURVEY / "epoch1" / f"{NAME}.eff")
        self.assertLessEqual(len(table), 30)          # OSSSSim n_bin_max
        mags = [m for m, _ in table]
        self.assertEqual(mags, sorted(mags))
        for m, e in table:
            self.assertAlmostEqual(e, P0 / 2 * math.erfc((m - M0) / W), places=4)
        self.assertAlmostEqual(table[0][1], P0, places=3)
        self.assertLess(table[-1][1], 0.01)            # below OSSSSim's eta floor

    def test_sample_a_depth(self):
        """p_det^3 is 40% at m = 28.6, the paper's Sample A limit (Section 10)."""
        table = lookup_table(SURVEY / "epoch1" / f"{NAME}.eff")
        mags = [m for m, _ in table]
        effs = [e for _, e in table]
        k = mags.index(28.6)
        self.assertAlmostEqual(effs[k] ** 3, 0.40, delta=0.01)

    def test_mag_error_curve(self):
        eff = read_keyed(SURVEY / "epoch1" / f"{NAME}.eff")
        e1, e2, e3, e4, _, _ = (float(x) for x in eff["mag_error"])
        alpha = math.log10(e2 / e1) / (e3 - 21.0)

        def sigma(m):
            if m <= e3:
                return e1 * 10 ** (alpha * (m - 21.0))
            return max(e2 - (m - e3) * e4, 0.0)

        # Background limited, through the Table 1 single-epoch uncertainties.
        self.assertAlmostEqual(alpha, 0.4, places=3)
        self.assertAlmostEqual(sigma(28.6), 0.11, delta=0.005)
        self.assertAlmostEqual(sigma(29.3), 0.21, delta=0.01)
        self.assertAlmostEqual(sigma(27.5), 0.04, delta=0.005)
        # The catalogue magnitude averages the three epochs.
        self.assertEqual([float(x) for x in eff["phot_frac"]], [0.0, 0.0, 1.0])

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

    def test_detections_are_sample_a(self):
        with (SURVEY / f"{NAME}.csv").open() as fh:
            rows = list(csv.DictReader(fh))
        by_comp = {}
        for row in rows:
            self.assertEqual(row["sample"], "A")
            by_comp.setdefault(row["comp"], set()).add(row["name"])
        self.assertEqual(by_comp, {"cold": COLD, "hot": HOT, "centaur": CENTAUR})
        # Eduardo et al. Section 7.1 dynamical classes.
        for row in rows:
            d, inc = float(row["d_bary"]), float(row["i"])
            if d < 30.0:
                self.assertEqual(row["comp"], "centaur", row["name"])
            elif 38.0 < d < 55.0 and inc <= 5.0:
                self.assertEqual(row["comp"], "cold", row["name"])
            else:
                self.assertEqual(row["comp"], "hot", row["name"])
        # Faintest Sample A member, and none of JPB02, JPB03 or 2015 GK56 (JPB23).
        self.assertEqual(max(float(r["m_f150w2"]) for r in rows), 28.60)
        self.assertFalse({"JPB02", "JPB03", "JPB23"} & {r["name"] for r in rows})


class SingleEpochSurveyLayoutTest(unittest.TestCase):
    def test_single_epoch_survey_roots(self):
        for name in SINGLE_EPOCH_SURVEYS:
            with self.subTest(survey=name):
                d = CHAR / name
                rows = pointing_rows(d / "pointings.list")
                self.assertTrue(rows)
                self.assertFalse((d / "survey.conf").exists())
                effs = {row[-1] for row in rows if row[-1].endswith(".eff")}
                self.assertTrue(effs)
                for eff in effs:
                    self.assertTrue((d / eff).is_file(), eff)

    def test_eff_files_match_pointings(self):
        for name in SINGLE_EPOCH_SURVEYS:
            with self.subTest(survey=name):
                d = CHAR / name
                referenced = {row[-1] for row in pointing_rows(d / "pointings.list")
                              if row[-1].endswith(".eff")}
                present = {p.name for p in d.glob("*.eff")}
                self.assertEqual(referenced, present)

    def test_detections_named_after_survey(self):
        for name in SINGLE_EPOCH_SURVEYS:
            with self.subTest(survey=name):
                d = CHAR / name
                self.assertEqual([p.name for p in d.glob("*.detections")],
                                 [f"{name}.detections"])
                objects = detection_objects(d / f"{name}.detections")
                self.assertTrue(objects)
                self.assertEqual(len(objects), len(set(objects)))
                cds = d / f"{name}.CDS"
                if cds.exists():
                    self.assertEqual(sorted(cds_objects(cds)), sorted(objects))

    def test_redetections_use_each_surveys_designation(self):
        pairs = pointing_rows(ROOT / "Duplicates")
        self.assertEqual(len(pairs), 6)
        cfeps = set(detection_objects(CHAR / "CFEPS_g" / "CFEPS_g.detections"))
        ossos = set(detection_objects(CHAR / "OSSOS_v11" / "OSSOS_v11.detections"))
        for ossos_name, cfeps_name in pairs:
            with self.subTest(ossos=ossos_name, cfeps=cfeps_name):
                self.assertIn(cfeps_name, cfeps)
                self.assertNotIn(ossos_name, cfeps)
                self.assertIn(ossos_name, ossos)
                self.assertNotIn(cfeps_name, ossos)

    def test_filters_have_colours(self):
        with (ROOT / "colour.toml").open("rb") as fh:
            colours = tomllib.load(fh)
        for name in SINGLE_EPOCH_SURVEYS:
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

    def test_single_epoch_surveys_load(self):
        from ossssim import SurveyCharacterization

        for name in SINGLE_EPOCH_SURVEYS:
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

        # A bright object on the field centre line of sight stays in the middle
        # strip, so it is found with probability (fill * p0)^3 = 0.67.
        draws = [flag(8.0) for _ in range(300)]
        found = [f for f, _ in draws if f >= 4]
        self.assertAlmostEqual(len(found) / len(draws), (FILL * P0) ** 3, delta=0.08)
        self.assertEqual({k for f, k in draws if f >= 4}, {f"{NAME}/{NAME}"})
        self.assertEqual(flag(16.0)[0], 0)


if __name__ == "__main__":
    unittest.main()
