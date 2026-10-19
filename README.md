# SurveyCharacterizations

Survey characterizations of trans-Neptunian object surveys, in the format
read by the OSSOS [SurveySimulator](https://github.com/ijiraq/SurveySimulator)
(`ossssim`). Per-detection biases are computed with
[survey_debias](https://github.com/ijiraq/survey_debias).

| Directory | Survey | Epochs | Filters |
|---|---|---|---|
| `CFEPS/` | CFEPS detections and blocks (OSSOS v11 release) | 1 | g, r, R |
| `OSSOS-MA/` | OSSOS plus the MA (mal, mah) blocks (v11) | 1 | r |
| `All_Surveys/` | OSSOS, CFEPS, HiLat, and MA: the full v11 sample | 1 | g, r, R |
| `All_r_Surveys/` | the r-band subset of `All_Surveys` (v11, no duplicates) | 1 | r |
| `Deep_Surveys/` | the 2015 OSSOS blocks (v11) | 1 | r |
| `OSSOS-B/` | OSSOS First Data Release, 2013–2014 blocks | 1 | r |
| `jwst_sample_a/` | JWST GO 1568 Sample A (Eduardo et al. 2026) | 3 | H (F150W2) |

All of them are under `characterization/`. Samples B and C of Eduardo et al.
will follow the same layout as `jwst_sample_b` and `jwst_sample_c`.

Orbit models (for example the OSSOS ModelUsed tables or the CFEPS L7 model)
are deliberately **not** kept here: they drive simulator runs and bias
calculations, but they are not characterizations. Keep them elsewhere and pass
them to `survey_debias` with `--model`.

## Layout

Every file that belongs to a survey lives in that survey's characterization
directory and is named after it:

```text
colour.toml                     colour table shared by every survey in a run
characterization/
  __README__.txt                OSSOS v11 release notes
  CFEPS/ OSSOS-MA/ All_Surveys/ All_r_Surveys/ Deep_Surveys/ OSSOS-B/
                                single-epoch surveys: pointings.list, .eff, detections
  jwst_sample_a/                multi-epoch survey root (no pointings.list of its own)
    jwst_sample_a.toml          survey_debias configuration
    jwst_sample_a.csv           detections (Eduardo et al. Sample A)
    survey.conf                 detections_required = 3
    epoch1/ epoch2/ epoch3/
      pointings.list            mosaic at the CADC shift-and-stack midpoint
      jwst_sample_a.eff         F150W2 detection efficiency (filter= H)
      JWST.csv                  JPL Horizons JWST barycentric vectors
```

- **Survey identity** is the directory name: detections are attributed to the
  key `<survey>/<block>`, for example `CFEPS/L3h-smooth`, `OSSOS-B/2013AE` or
  `jwst_sample_a/jwst_sample_a`.
- **Multi-epoch surveys** have no `pointings.list` in the survey directory.
  `ossssim` loads every child directory that holds one as an epoch, and
  `survey.conf` sets `detections_required`, the number of epochs at which a
  simulated object must be found.
- **`colour.toml` stays at the root** because SurveySimulator treats colour as
  a property of the simulated object, not of a survey: each simulated object
  carries one colour per filter letter, and each `.eff` file only picks the
  letter. All surveys in a run therefore share one colour table, and different
  instrument bands need different filter letters.

## OSSOS v11 and CFEPS

`CFEPS/`, `OSSOS-MA/`, `All_Surveys/`, `All_r_Surveys/` and `Deep_Surveys/`
are verbatim copies of the OSSOS v11 release characterizations (see
`characterization/__README__.txt`), kept under their original names so
detection keys match other OSSOS work. Do not edit them.

## OSSOS-B

`OSSOS-B/` is the OSSOS First Data Release characterization: the 2013–2014
blocks E, O, L and H with the detections released with Bannister et al.
(2016), AJ 152, 70,
[doi:10.3847/0004-6256/152/3/70](https://doi.org/10.3847/0004-6256/152/3/70).
It is kept, unmodified, so the figures of that paper can be reproduced; for
new work use the v11 sets, which supersede it. Its `.eff` files are identical
to the v11 ones; its pointings and detection lists are those of the first
release.

## JWST Sample A

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

## Installation

Both packages go into the same Python environment (≥ 3.12):

```bash
cd ../SurveySimulator && SDKROOT=$(xcrun --show-sdk-path) pip install .
cd ../survey_debias && pip install .
```

`SDKROOT` is only needed with an older Homebrew GCC whose built-in macOS SDK
path no longer exists. `survey_debias` must be a version that reads
`survey.conf` for multi-epoch surveys.

## Computing the bias

Run from the repository root; paths in the TOML are relative to `--root`:

```bash
survey_debias --survey characterization/jwst_sample_a/jwst_sample_a.toml --root . --dry-run
survey_debias --survey characterization/jwst_sample_a/jwst_sample_a.toml --root . \
    --model /path/to/OSSOS/models
```

`model_ae` surveys such as Sample A need an orbit model for p(a, e | r, i);
see the survey_debias README for the OSSOS model tables.

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
