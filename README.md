# SurveyCharacterizations

Survey characterizations of trans-Neptunian object surveys, in the format
read by the OSSOS [SurveySimulator](https://github.com/ijiraq/SurveySimulator)
(`ossssim`). Per-detection biases are computed with
[survey_debias](https://github.com/ijiraq/survey_debias).

| Directory | Survey | Epochs | Filters |
|---|---|---|---|
| `CFEPS_g/` | the g-band CFEPS blocks (v11 release) | 1 | g |
| `CFEPS_r/` | the CFEPS L3h block (v11 release) | 1 | r |
| `CFEPS_presurvey/` | the CFEPS presurvey block (v11 release) | 1 | R |
| `HiLat/` | the High Latitude survey (v11 release) | 1 | r |
| `MA/` | the MA (mal, mah) blocks (v11 release) | 1 | r |
| `OSSOS_v11/` | OSSOS, all blocks (v11 release) | 1 | r |
| `jwst_sample_a/` | JWST GO 1568 Sample A (Eduardo et al. 2026) | 3 | H (F150W2) |
| `N26/` | HST GO-9433 heliostack (Napier et al. 2026) | 1 | t (F606W STMAG) |

All of them are under `characterization/`. Sample B of Eduardo et al. (two
epochs at S/N > 8, a different efficiency) would follow the same layout as
`jwst_sample_b`. The paper has no sample C: the "C" in its Table 1 marks
cold classicals.

Orbit models (for example the OSSOS ModelUsed tables or the CFEPS L7 model)
are deliberately **not** kept here: they drive simulator runs and bias
calculations, but they are not characterizations. Keep them elsewhere and pass
them to `survey_debias` with `--model`.

## Layout

Every file that belongs to a survey lives in that survey's characterization
directory and is named after it:

```text
colour.toml                     colour table shared by every survey in a run
Duplicates                      OSSOS and CFEPS designations of the 6 objects in both surveys
characterization/
  __README__.txt                OSSOS v11 release notes
  CFEPS_g/ CFEPS_r/ CFEPS_presurvey/ HiLat/ MA/ OSSOS_v11/
                                single-epoch surveys: pointings.list, .eff,
                                <survey>.detections, optional <survey>.CDS
  jwst_sample_a/                multi-epoch survey root (no pointings.list of its own)
    README.md                   where every value came from, and the checks
    jwst_sample_a.toml          survey_debias configuration
    jwst_sample_a.csv           detections (Eduardo et al. Sample A)
    survey.conf                 detections_required = 3
    footprints.csv              MAST footprints of the F150W2 SW exposures
    make_pointings.py           writes the epoch pointings from footprints.csv
    epoch1/ epoch2/ epoch3/
      pointings.list            outline of the three NIRCam strips, fill 0.91
      jwst_sample_a.eff         F150W2 single-epoch efficiency (filter= H)
      JWST.csv                  JPL Horizons JWST barycentric vectors
  N26/                          single-epoch survey with its own observer and config
    README.n26                  where every value came from
    N26.toml                    survey_debias configuration
    N26.csv                     detections (Napier et al. Tables 1 and 2)
    pointings.list              ACS/WFC mosaic at the 15-day stack midpoint
    N26.eff                     F606W STMAG detection efficiency (filter= t)
    HST.csv, HST.jpl_batch      JPL Horizons HST barycentric vectors and request
```

- **Survey identity** is the directory name: detections are attributed to the
  key `<survey>/<block>`, for example `CFEPS_r/L3h-smooth`, `OSSOS_v11/2013AE`
  or `jwst_sample_a/jwst_sample_a`.
- **One survey per directory**: every `.eff` named in `pointings.list` is in
  the same directory, and every `.eff` in the directory is named there.
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

`CFEPS_g/`, `CFEPS_r/`, `CFEPS_presurvey/`, `HiLat/`, `MA/` and `OSSOS_v11/`
split the OSSOS v11 release characterizations (see
`characterization/__README__.txt`) into one directory per survey; CFEPS is
split by filter. The `.eff` files, pointings and detection rows are those of
the release, except that the six CFEPS objects re-detected by OSSOS carry
their CFEPS designations in `CFEPS_g/` (the release used the OSSOS ones).
Detection keys use the new directory names (for example `CFEPS_g/L3f-smooth`
rather than the release's `CFEPS/L3f-smooth`).

The six re-detected objects are in both `CFEPS_g/` (CFEPS designation and
circumstances) and `OSSOS_v11/` (OSSOS designation and circumstances);
`Duplicates` maps one name to the other. A run that combines the two surveys
counts them twice.

## JWST Sample A

| | |
|---|---|
| Field | the three NIRCam strips of the 20-tile mosaic (178.9 arcmin², one polygon per epoch), fill 0.91 |
| Epochs (JD) | 2459969.32121, 2459973.96788, 2459979.90857 |
| Detection | found at all 3 epochs (`detections_required = 3`) |
| Rates | 0.03–2.8 ″/hr (epoch 1), 0.03–1.4 ″/hr (epochs 2, 3), any direction |
| Efficiency per epoch | 0.96/2 erfc((m − 28.92)/0.61), F150W2, as a lookup table |
| Sample A depth | p_det³ = 40% at m = 28.6 |
| Magnitude limit | 30.5 (F150W2) |
| Photometric error | 0.11 × 10^(0.4(m − 28.6)) per epoch; catalogue magnitude is the 3-epoch mean |
| Detections | 20 (12 cold, 7 hot, 1 centaur) |

F150W2 is labelled with the filter letter `H` (close to 2MASS H).
`colour.toml` sets `"H-g" = -1.7`, so r − H = 1.0, i.e. m_r = m_F150W2 + 1
as in the paper.

`characterization/jwst_sample_a/README.md` gives the source of every value and
describes the checks against Eduardo et al. (2026) that set them. In summary:

- The sample is right.
- The efficiency parameters are the per-epoch ones, cubed by the 3-epoch
  requirement, but the old `.eff` used a tanh where the paper has an erfc.
- The old 0.05 deg² square overstated the three-epoch area by 40–70%.
- The photometric errors were not the paper's.
- JPB08 is cold, and JPB09 is a centaur.

## N26 (HST heliostack)

| | |
|---|---|
| Field | 14:07:53.33 −11:21:38, 400″ × 600″ ACS/WFC mosaic |
| Epoch (JD) | 2452672.8585, midpoint of the 15-day heliostack |
| Detection | single stack |
| Rates | 0.05–6.4 ″/hr, any direction |
| Efficiency | 1.0/2 (1 − tanh((m − 29.21)/0.22)), F606W STMAG |
| Magnitude limit | 29.8 (F606W STMAG) |

F606W STMAG is labelled with the filter letter `t`: it is close to r but is
neither r nor w. `colour.toml` sets `"t-g" = "r-g" + 0.3`, i.e.
r_AB = STMAG − 0.3. Provenance of every number is in
`characterization/N26/README.n26`.

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
