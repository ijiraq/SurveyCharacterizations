"""Checks on the MPC complete-survey characterization (characterization/MPC)."""
from __future__ import annotations

import csv
import math
import re
import sys
import tomllib
import unittest
from pathlib import Path

from test_characterization import poly_pointings, read_keyed

ROOT = Path(__file__).resolve().parents[1]
NAME = "MPC"
SURVEY = ROOT / "characterization" / NAME
sys.path.insert(0, str(SURVEY))

import make_mpc  # noqa: E402
import make_pointings  # noqa: E402


def load_spec() -> dict:
    with (SURVEY / f"{NAME}.toml").open("rb") as fh:
        return tomllib.load(fh)


def sample_rows() -> list[dict]:
    lines = [x for x in (SURVEY / f"{NAME}.csv").read_text().splitlines()
             if x.strip() and not x.startswith("#")]
    return list(csv.DictReader(lines))


class MPCLayoutTest(unittest.TestCase):
    def test_single_epoch_tiles(self):
        self.assertFalse((SURVEY / "survey.conf").exists())
        tiles = poly_pointings(SURVEY / "pointings.list")
        self.assertGreater(len(tiles), 0)
        self.assertLessEqual(len(tiles), 2000)
        for head, verts in tiles:
            self.assertEqual(head[0], "poly")
            self.assertEqual(int(head[1]), len(verts))
            self.assertLess(len(verts), 41)
            self.assertEqual(float(head[5]), 1.0)
            self.assertEqual(head[6], "500")
            self.assertRegex(head[7], rf"^{NAME}_\d{{3}}\.eff$")
        self.assertEqual({p.name for p in SURVEY.glob("*.eff")}, {h[7] for h, _ in tiles})

    def test_tiles_cover_every_ecliptic_longitude(self):
        tiles = poly_pointings(SURVEY / "pointings.list")
        jds = sorted({float(h[4]) for h, _ in tiles})
        self.assertEqual(len(jds), 36)
        self.assertLess(jds[-1] - jds[0], 366.0)

    def test_efficiency_imitates_h_limit(self):
        spec = load_spec()["survey"]
        for path in sorted(SURVEY.glob("*.eff")):
            with self.subTest(eff=path.name):
                eff = read_keyed(path)
                self.assertEqual(eff["filter"], [spec["filter"]])
                self.assertEqual([float(x) for x in eff["mag_error"]], [0.001, 0, 99, 0, 99, 0])
                self.assertEqual([float(x) for x in eff["track_frac"]], [1.0, 99.0, 0.0])
                lines = path.read_text().splitlines()
                earth = next(x for x in lines if x.startswith("# earth:")).split()
                vals = dict(zip(earth[2::2], (float(x) for x in earth[3::2])))
                ranges = []
                for k, line in enumerate(lines):
                    if line.startswith("rates="):
                        lo, hi = (float(x) for x in line.split("=")[1].split())
                        _, m0, _, alpha = (float(x) for x in lines[k + 2].split("=")[1].split())
                        ranges.append((lo, hi, m0, alpha))
                self.assertLessEqual(len(ranges), 10)
                self.assertEqual(ranges[0][0], 0.0)
                for (_, hi, _, _), (lo, _, _, _) in zip(ranges[:-1], ranges[1:]):
                    self.assertAlmostEqual(hi, lo, places=4)
                self.assertLess(abs(vals["r_e="] - 1.0), 0.02)
                self.assertGreater(vals["phase="], 0.0)
                for r in (38.0, 41.0, 44.0, 47.0, 51.0, 55.0):
                    rate = float(make_pointings.opposition_rate(r, vals["r_e="], vals["v_e="]))
                    lo, hi, m0, alpha = next(x for x in ranges if x[0] <= rate <= x[1])
                    want = (spec["h_lim"] + 5.0 * math.log10(r * (r - vals["r_e="]))
                            + vals["phase="])
                    self.assertAlmostEqual(m0 + alpha * rate, want, delta=0.01)

    def test_filter_colours(self):
        with (ROOT / "colour.toml").open("rb") as fh:
            colours = tomllib.load(fh)
        cold, default = colours["cold"], colours["default"]
        self.assertEqual(load_spec()["survey"]["filter"], "r")
        for group in (cold, default):
            gr = -group["r-g"]
            self.assertAlmostEqual(group["V-g"], -0.59 * gr - 0.01, delta=0.01)

    def test_h_is_r_band(self):
        rows = sample_rows()
        meta = [x for x in (SURVEY / f"{NAME}.csv").read_text().splitlines() if x.startswith("#")]
        m = re.search(r"h_mpc - ([0-9.]+)", next(x for x in meta if x.startswith("# h:")))
        offset = float(m.group(1))
        lines = [x for x in (SURVEY / "MPC_photometry.csv").read_text().splitlines()
                 if not x.startswith("#")]
        measured = {r["desig"]: float(r["H_r"]) for r in csv.DictReader(lines)}
        sources = {"photometry", "mpc_offset", "mpc_offset_photometry_rejected"}
        for r in rows:
            self.assertIn(r["h_source"], sources)
            if r["h_source"] == "photometry":
                self.assertAlmostEqual(float(r["h"]), measured[r["desig"]], places=3)
            else:
                self.assertAlmostEqual(float(r["h"]), float(r["h_mpc"]) - offset, places=3)
        spec = load_spec()["survey"]
        bright = [r for r in rows if r["comp"] == "cold" and float(r["h"]) < spec["h_lim"]]
        self.assertTrue(all(r["h_source"] == "photometry" for r in bright))

    def test_spec_paths_in_survey_directory(self):
        spec = load_spec()["survey"]
        self.assertEqual(spec["surveys"], [NAME])
        self.assertEqual(spec["format"], "complete")
        for key in ("detections_relpath", "duplicates_relpath"):
            self.assertEqual(Path(spec[key]).parent, SURVEY.relative_to(ROOT), key)

    def test_sample_columns_and_cut(self):
        spec = load_spec()
        rows = sample_rows()
        self.assertGreater(len(rows), 0)
        for col in spec["detection_columns"].values():
            self.assertIn(col, rows[0])
        names = [r["desig"] for r in rows]
        self.assertEqual(len(names), len(set(names)))
        cut = spec["cut"]
        cold = [r for r in rows if r["comp"] == "cold"]
        for r in cold:
            self.assertEqual((r["cl"], r["p"]), ("cla", "m"))
            self.assertLess(float(r["ifree"]), cut["ifree_max"])
            self.assertIn(r["ifree_source"], ("huang2022_table", "huang2022_grid"))
            self.assertGreaterEqual(float(r["a"]), cut["a_min"])
            self.assertLessEqual(float(r["a"]), cut["a_max"])
        bright = [r for r in cold if float(r["h"]) < spec["survey"]["h_lim"]]
        self.assertGreater(len(bright), 0)
        self.assertTrue(all(r["secure"] == "1" for r in bright))

    def test_vetting_applied(self):
        vetting = make_mpc.read_vetting(SURVEY / "MPC.vetting")
        self.assertGreater(len(vetting), 0)
        by_desig = {r["desig"]: r for r in sample_rows()}
        huang = {"0": "cla", "-1": "sca"}
        for desig, (cl, _) in vetting.items():
            row = by_desig[desig]
            want = cl if row["huang_reso"] == "" else huang.get(row["huang_reso"], "res")
            self.assertEqual(row["cl"], want, desig)

    def test_duplicates_name_real_objects(self):
        names = {r["desig"] for r in sample_rows()}
        release = {obj for _, obj, _ in make_mpc.release_names(ROOT / "characterization")}
        for raw in (SURVEY / f"{NAME}.duplicates").read_text().splitlines():
            words = raw.split("#", 1)[0].split()
            if len(words) >= 2:
                self.assertIn(words[0], names)
                self.assertIn(words[1], release)

    def test_unpack_designations(self):
        self.assertEqual(make_mpc.unpack_provisional("K13SA0A"), "2013SA100")
        self.assertEqual(make_mpc.unpack_provisional("K14UM5E"), "2014UE225")
        self.assertEqual(make_mpc.unpack_provisional("J98W00Z"), "1998WZ")
        self.assertEqual(make_mpc.unpack_provisional("C9878"), "129878")
        self.assertEqual(make_mpc.unpack_epoch("K2669"), "2026-06-09")


def _ossssim():
    try:
        import ossssim
    except ImportError:
        return None
    return ossssim if hasattr(ossssim, "read_survey_conf") else None


@unittest.skipIf(_ossssim() is None, "needs ossssim with survey roots")
class MPCSimulatorTest(unittest.TestCase):
    def test_loads_as_single_epoch_survey(self):
        from ossssim import SurveyCharacterization

        s = SurveyCharacterization.from_directory(str(SURVEY))
        self.assertEqual(s.n_epochs, 1)
        self.assertEqual(s.detections_required, 1)
        self.assertEqual(len(s), len(poly_pointings(SURVEY / "pointings.list")))
        self.assertEqual({p.survey for p in s}, {NAME})


if __name__ == "__main__":
    unittest.main()
