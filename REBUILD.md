# Rebuilding every bias product from scratch

Use this after changing SurveySimulator (`ossssim`), `survey_debias`, or any
characterization. Commands run from the repository root, so `--root .`
(the TOML comments in `joint/` use `--root ..` because they run from `joint/`).

There are two separate products:

| Product | Tool | Surveys | Files (all gitignored) |
|---|---|---|---|
| Per-detection bias | `survey_debias` | surveys with a `<survey>.toml`: `jwst_sample_a`, `N26`, `DEEP` | `bias_grid*.csv` (cache), `<survey>.bias.csv`, `<survey>.detections-full`, `check_plots/` |
| Population efficiency eps(H) for the joint fits | `survey_debias_eff` | `OSSOS_v11`, `CFEPS_r`, `CFEPS_g`, `CFEPS_presurvey`, `MA`, `jwst_sample_a`, `N26`, `MPC` (analytic) | `plants_<tag>.npz`, `eps_<tag>.csv`, `eps_<tag>_x_<other>.csv` |

Both eps tags are used: `cold` (ecliptic i < 5°, `joint/cold_*.toml`) and
`coldfree` (free inclination < 4°, `joint/coldfree_*.toml`).

## 0. Check the installs

```bash
python -c "import ossssim, survey_debias; print(ossssim.__file__, survey_debias.__file__)"
(cd ../SurveySimulator && python -m pytest -q tests/test_observer_frame.py)   # equinox fix present
python -m unittest discover -s tests
```

The orbit model is not in the repository (`Models/` is gitignored):

```bash
mkdir -p Models/OSSOS
for c in Classical Detached Inner Plutinos Scattering Twotinos; do
  curl -fsSL -o Models/OSSOS/$c-ModelUsed.dat \
    https://raw.githubusercontent.com/OSSOS/OSSOS_Models/$c/$c/ModelUsed-check-8.66.dat
done
M=Models/OSSOS/Classical-ModelUsed.dat
```

## 1. Remove every cached result

`bias_grid*.csv` caches are reused automatically, so stale cells survive a re-run
unless they are deleted.

```bash
rm -f characterization/*/bias_grid*.csv characterization/*/plants_*.npz \
      characterization/*/eps_*.csv characterization/*/*.bias.csv \
      characterization/*/*.detections-full
rm -rf check_plots joint/joint_*/
```

## 2. Re-check the excluded blocks

The joint fits leave out OSSOS 2015BS and CFEPS L3s-smooth and L3f-smooth because
their characterizations recovered < 30% of their own detections. Detos1
rejected the detections' direction of motion, which is consistent with the
observer-frame bug fixed in ossssim on 2026-10-20. Check them again:

```bash
survey_debias_eff --root . --model $M --check-detections \
    --survey OSSOS_v11 --survey CFEPS_g --survey CFEPS_r --survey CFEPS_presurvey --survey MA
```

Lines marked `<-- low` still fail. If 2015BS, L3s-smooth or L3f-smooth no
longer fail, drop their `--exclude-block` below and their `exclude_blocks`
in `joint/*.toml`.

```bash
X="--exclude-block 2015BS --exclude-block L3s-smooth --exclude-block L3f-smooth"
```

## 3. Per-detection biases

```bash
for s in jwst_sample_a N26 DEEP; do
  survey_debias --survey characterization/$s/$s.toml --root . --dry-run
done
survey_debias --survey characterization/jwst_sample_a/jwst_sample_a.toml --root . --model Models/OSSOS
survey_debias --survey characterization/N26/N26.toml --root . --model Models/OSSOS
survey_debias --survey characterization/DEEP/DEEP.toml --root .
```

The default is `--target 5000` detections per cell. DEEP has 110 cells over 29
epochs: `--target 100` took 6 minutes, so expect several hours at the default.
Use `--target 200` to check the setup first. Each run prints `sum 1/bias`.

## 4. Population efficiency tables

Order matters: `--overlap-with SURVEY` needs that survey's plants, so each
survey comes after the ones it overlaps. The MPC is analytic and goes last.

### Tag `cold` (ecliptic i < 5°)

```bash
C="--root . --model $M --i-max 5"
survey_debias_eff $C --survey OSSOS_v11 $X
survey_debias_eff $C --survey CFEPS_r $X --overlap-with OSSOS_v11
survey_debias_eff $C --survey CFEPS_g $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r
survey_debias_eff $C --survey CFEPS_presurvey $X --overlap-with OSSOS_v11 \
    --overlap-with CFEPS_r --overlap-with CFEPS_g
survey_debias_eff $C --survey MA $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r \
    --overlap-with CFEPS_g --overlap-with CFEPS_presurvey
survey_debias_eff --root . --model $M --i-max 5.0001 --r-min 38 --r-max 55 \
    --survey jwst_sample_a $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r --overlap-with CFEPS_g
survey_debias_eff $C --survey N26 $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r \
    --overlap-with CFEPS_g --overlap-with CFEPS_presurvey --overlap-with MA \
    --overlap-with jwst_sample_a
survey_debias_eff --root . --complete characterization/MPC/MPC.toml --model $M --i-max 5
```

### Tag `coldfree` (i_free < 4°; plants keep the model's free inclination)

```bash
C="--root . --model $M --tag coldfree --ifree-max 4"
survey_debias_eff $C --survey OSSOS_v11 $X
survey_debias_eff $C --survey CFEPS_r $X --overlap-with OSSOS_v11
survey_debias_eff $C --survey CFEPS_g $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r
survey_debias_eff $C --survey CFEPS_presurvey $X --overlap-with OSSOS_v11 \
    --overlap-with CFEPS_r --overlap-with CFEPS_g
survey_debias_eff $C --survey MA $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r \
    --overlap-with CFEPS_g --overlap-with CFEPS_presurvey
survey_debias_eff $C --survey jwst_sample_a $X --r-min 38 --r-max 55 \
    --overlap-with OSSOS_v11 --overlap-with CFEPS_r --overlap-with CFEPS_g
survey_debias_eff $C --survey N26 $X --overlap-with OSSOS_v11 --overlap-with CFEPS_r \
    --overlap-with CFEPS_g --overlap-with CFEPS_presurvey --overlap-with MA \
    --overlap-with jwst_sample_a
survey_debias_eff --root . --complete characterization/MPC/MPC.toml --model $M \
    --tag coldfree --ifree-max 4
```

`--plants` defaults to 20000 aimed plants per block. Without `--reuse` every
run plants again.

### DEEP (not yet in any joint fit)

DEEP's eps can be built the same way. Its cold sample is
42.4 ≤ a ≤ 47.7 au, i < 5°, q ≥ 37 au (`DEEP.csv` `comp = cold`).
`survey_debias_eff` has no q cut, so this is the closest match:

```bash
survey_debias_eff --root . --model $M --i-max 5 --a-min 42.4 --a-max 47.7 \
    --survey DEEP --overlap-with OSSOS_v11 --overlap-with CFEPS_r --overlap-with CFEPS_g \
    --overlap-with CFEPS_presurvey --overlap-with MA --overlap-with jwst_sample_a --overlap-with N26
```

Using it in a fit needs a `[[survey]]` entry and a `priority` position in the
`joint/*.toml` files. That has not been added; DEEP's efficiency is still
provisional (`characterization/DEEP/README.deep`).

## 5. Joint fits

```bash
cd joint
for cfg in cold_jwst_ossos cold_jwst_ossos_cfeps cold_all cold_all_mpc coldfree_all coldfree_all_mpc; do
  survey_debias_joint --config $cfg.toml          # writes joint/joint_<cfg>/
done
```

## 6. Checks

- `python -m unittest discover -s tests` passes.
- No block prints `<-- low` in step 2, apart from those deliberately excluded.
- Compare each `joint_*/` result with the previous one. Expect shifts from both
  2026-10-20 fixes:
  - the equinox observer frame, which affects code-500 blocks observed within a month of
    an equinox;
  - the single-attempt aiming, which affects cells whose inclination is close to the field's
    ecliptic latitude, in every survey.
