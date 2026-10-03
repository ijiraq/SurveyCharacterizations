# JWST / Roman TNO follow-up

Standalone project that **uses** [SurveySimulator](https://github.com/ijiraq/SurveySimulator)
(`ossssim`). Characterization, Sample A grid-bias wrappers, and the GBTDS
sampling catalog live here so they are not distributed with the simulator.

Clone this repository next to `SurveySimulator`, or set `SSIM_ROOT` to that
checkout. The CFEPS L7 model is **not** vendored; the catalog generator reads
`F95/tests/Models/L7model-3.0-9.0` from SSim.

```text
jwst-tno-followup/          this repo
SurveySimulator/            sibling clone (or $SSIM_ROOT)
```

## GBTDS sampling catalog

CFEPS L7 sky positions around the Roman GBTDS pointing for the JWST follow-up
proposal. Detectability is not applied.

| | |
|---|---|
| Pointing | `13:52:25.52` `−11:01:25.3` (ICRS) |
| Epoch | 2027-05-01 (geocenter, barycentric ICRF) |
| Radius | 10° |
| Model | `$SSIM_ROOT/F95/tests/Models/L7model-3.0-9.0` |
| In field | 3676 / 106074 (median `helio_dist` 44.95 AU) |

```bash
pip install -r requirements.txt
python scripts/gbtds_field_positions.py
python -m unittest discover -s tests
```

- `data/gbtds_radec_helio_2027may01.csv` — `ra`, `dec`, `helio_dist` for sampling
- `data/gbtds_radec_helio_2027may01.ecsv` — same rows plus `delta`, `sep_deg`, and L7 elements

Propagation matches Detos1 (`sqrt(gmb)` mean motion, `pos_cart` + `RADECeclXV`)
without loading a survey characterization.

## Sample A grid bias

`scripts/grid_bias.py` and `scripts/compute_grid_bias.py` wrap
`ossssim.grid_bias` for the JWST Sample A mosaic. The Horvitz–Thompson helpers
remain in SurveySimulator; this repo only holds the survey characterization and
the thin wrappers.

Grid-bias runs need an installed `ossssim` (Fortran extension) plus
`characterization/` in this tree:

```bash
export SSIM_ROOT=/path/to/SurveySimulator
python scripts/compute_grid_bias.py
```

## Licence

Released under the European Union Public Licence (EUPL), same as
SurveySimulator. See `LICENCE.txt`.
