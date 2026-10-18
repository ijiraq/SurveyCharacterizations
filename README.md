# JWST TNO characterizations (Eduardo et al. 2026)

Survey characterizations of the JWST GO 1568 trans-Neptunian object samples
of Eduardo et al. 2026 (AJ 172, 187), in the format read by the OSSOS
[SurveySimulator](https://github.com/ijiraq/SurveySimulator) (`ossssim`).
Per-detection biases are computed with
[survey_debias](https://github.com/ijiraq/survey_debias).

Sample A is characterized. Samples B and C will follow the same layout as
`jwst_sample_b` and `jwst_sample_c`.

The ground-based CFEPS/OSSOS characterizations from SurveySimulator-Data are
included alongside them (see
[Characterizations from SurveySimulator-Data](#characterizations-from-surveysimulator-data)),
so JWST and ground-based surveys can be simulated from one tree.

## Layout

Every file that belongs to a survey lives in that survey's characterization
directory and is named after it:

```text
colour.toml                     colour table shared by every survey in a run
characterization/
  jwst_sample_a/                survey root (no pointings.list of its own)
    jwst_sample_a.toml          survey_debias configuration
    jwst_sample_a.csv           detections (Eduardo et al. Sample A)
    survey.conf                 detections_required = 3
    epoch1/ epoch2/ epoch3/
      pointings.list            mosaic at the CADC shift-and-stack midpoint
      jwst_sample_a.eff         F150W2 detection efficiency (filter= H)
      JWST.csv                  JPL Horizons JWST barycentric vectors
```

`colour.toml` stays at the root because SurveySimulator treats colour as a
property of the simulated object, not of a survey: each simulated object
carries one colour per filter letter, and each `.eff` file only picks the
letter. All surveys in a run therefore share one colour table, and different
instrument bands need different filter letters.

`characterization/jwst_sample_a` is a multi-epoch survey root: `ossssim`
loads every child directory that holds a `pointings.list` as one epoch, and
`survey.conf` makes a simulated object count as detected only if it is found
at all three epochs. Detections are attributed to the key
`jwst_sample_a/jwst_sample_a` (survey/block).

## Sample A

| | |
|---|---|
| Field | 13:57:33 −10:51:55 (ICRS), 0.05 deg² square mosaic |
| Epochs (JD) | 2459969.32118, 2459973.96785, 2459979.90854 |
| Detection | found at all 3 epochs (`detections_required = 3`) |
| Rates | 0.03–8.66 ″/hr, any direction |
| Efficiency | 0.96/2 (1 − tanh((m − 28.92)/0.61)), F150W2 |
| Magnitude limit | 30.5 (F150W2) |
| Photometric error | 0.05 at m = 28.2 rising to 0.5 at m = 29.5 |

F150W2 is labelled with the filter letter `H` (close to 2MASS H).
`colour.toml` sets `"H-g" = -1.7`, so r − H = 1.0, i.e. m_r = m_F150W2 + 1
as in the paper.

## Characterizations from SurveySimulator-Data

These directories are verbatim copies of
`SurveySimulator-Data/Characterizations/` (the OSSOS v11 release data; see
`characterization/__README__.txt`). They are single-epoch survey roots, each
with its own `pointings.list`, `.eff` blocks, and detections files, and are
kept under their original names so detection keys (for example
`CFEPS/L3h-smooth`) match other OSSOS work. Do not edit them here.

| Directory | Contents | Filters |
|---|---|---|
| `CFEPS/` | CFEPS detections and blocks | g, r, R |
| `OSSOS-MA/` | OSSOS plus the MA (mal, mah) blocks | r |
| `All_Surveys/` | OSSOS, CFEPS, HiLat, and MA: the full sample | g, r, R |
| `All_r_Surveys/` | the r-band subset of `All_Surveys` (no duplicates) | r |
| `Deep_Surveys/` | the 2015 OSSOS blocks | r |

Their filters are covered by the shared `colour.toml`. A `survey_debias`
configuration selects which surveys it uses with `surveys = [...]`, so these
directories do not affect a Sample A run.

## Installation

Both packages go into the same Python environment (≥ 3.12):

```bash
# SurveySimulator with multi-epoch survey roots
cd ../SurveySimulator && git switch cursor/multi-epoch-survey-root-5dc8
SDKROOT=$(xcrun --show-sdk-path) pip install .     # SDKROOT: macOS + Homebrew gfortran
# survey_debias
cd ../survey_debias && pip install .
```

`model_ae` biases need the OSSOS model tables in `Models/OSSOS/`
(gitignored); see the survey_debias README for the download.

## Computing the bias

Run from the repository root; paths in the TOML are relative to `--root`:

```bash
survey_debias --survey characterization/jwst_sample_a/jwst_sample_a.toml --root . --dry-run
survey_debias --survey characterization/jwst_sample_a/jwst_sample_a.toml --root . \
    --model Models/OSSOS
```

Outputs, in `characterization/jwst_sample_a/`:

- `jwst_sample_a.bias.csv`: the bias of each detection, with its Monte Carlo
  error.
- `jwst_sample_a.detections-full`: CFEPS-format detections with their bias.
- `bias_grid_rih.csv`: per-cell cache, reused on re-runs (gitignored; delete
  it to recompute).

Check plots go to `check_plots/` at the root (gitignored).

The default `--target` is 5000 detections per cell; `--target 200` takes
about 20 s for Sample A and is enough to check the setup.

## Tests

```bash
python -m unittest discover -s tests
```

The simulator tests are skipped when the installed `ossssim` cannot load
multi-epoch survey roots.

## Licence

Released under the European Union Public Licence (EUPL), same as
SurveySimulator. See `LICENCE.txt`.
