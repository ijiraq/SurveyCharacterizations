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
| `MPC/` | the MPC orbit catalogue, complete for cold classicals with H_r < 5.42 | 1 (36 tiles) | r |
| `DEEP/` | DEEP B1 linked sample, 110 TNOs (Smotherman et al. 2024); efficiency provisional | 29 (any 4) | r (DECam VR on PS1 r) |

All of them are under `characterization/`. 

## Layout

Every file that belongs to a survey lives in that survey's characterization
directory and is named after it:

```text
colour.toml                     TNO average colours table shared by every survey
Duplicates                      Objects that appear in multiple surveys. 
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
  MPC/                          complete survey (all cold classicals below an H limit)
    README.mpc                  completeness limit, H scale, vetting, ossssim validation
    MPC.toml                    complete-survey spec (filter r, h_lim, sigma_h, cut)
    MPC.csv                     MPCORB orbits with 30 <= a <= 150 au, MPC H <= 7; h = H_r
    mpc_photometry.py, obs/     MPC observations -> measured H_r (MPC_photometry.csv)
    MPC.duplicates, MPC.vetting cross-IDs with OSSOS/CFEPS; resonance vetting
    make_mpc.py                 MPCORB.DAT (not committed) -> MPC.csv, completeness, calibration
    vet_resonances.py           10 Myr REBOUND integrations behind MPC.vetting
    make_pointings.py           pointings.list and MPC_<lon>.eff (ossssim imitation)
    validate_ossssim.py         checks the tiles reproduce H_r < h_lim
  DEEP/                         multi-epoch survey: 29 long stares, detection on any 4
    README.deep                 sources, checks, approximations, what is still needed
    DEEP.toml                   survey_debias configuration
    DEEP.csv                    detections (DEEP VI Table 2)
    survey.conf                 detections_required = 4
    make_deep.py                source/ -> DEEP.csv and epoch*/
    source/                     VizieR table, DEEP VI Table 1, NOIRLab archive stares
    epoch1/ ... epoch29/
      pointings.list            DECam outline at the stare's boresight, fill 0.860
      DEEP.eff                  single-night efficiency (filter= r), provisional
```

- **Survey identity** is the directory name: detections are attributed to the
  key `<survey>/<block>`, for example `CFEPS_r/L3h-smooth`, `OSSOS_v11/2013AE`
  or `jwst_sample_a/jwst_sample_a`.
- **One survey per directory**: every `.eff` named in `pointings.list` is in
  the same directory, and every `.eff` in the directory is named there.
- **Multi-epoch surveys** have no `pointings.list` in the survey directory.
  `ossssim` loads every child directory that holds a pointing.list file as an epoch, and
  `survey.conf` sets `detections_required`, the number of epochs at which a
  simulated object must be found.  The fortran branch also supports multiple epochs.
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
as documented in Eduardo et al. (2026).

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

## MPC (complete survey)

| | |
|---|---|
| Sample | MPCORB snapshot 2026-10-19; cold classicals: 42.4 ≤ a ≤ 47.7 au, free inclination i_free < 4° (Huang et al. 2022), q ≥ 37 au, non-resonant (Huang et al.'s classification, else 10 Myr integrations) |
| Completeness | H_r < 5.42 (MPC H < 5.5): no cold discovery brighter than that since 2014 (the whole Kuiper belt is *not* complete there) |
| Detections | 5, all with H_r measured from their MPC photometry |
| H scale | r (PS1/Rubin): `MPC.csv` `h` is H_r, `h_mpc` the MPC's own H |
| survey_debias | `format = "complete"`: eps = 0.714 × 1(H_r < 5.42), analytic (`survey_debias_eff --complete --ifree-max 4`) |
| ossssim | 36 ecliptic tiles (\|β\| ≤ 8°) at opposition, per-tile `.eff` with a rate-dependent limit; recovers H_r < 5.42 to 0.005 mag (edge 10–90% in 0.15 mag) |

MPCORB's H is the MPC's own, from mixed bands with its own corrections, so
`MPC.csv` carries H_r instead. `mpc_photometry.py` downloads each candidate's
observations from the MPC, reduces them to H in each reported band from the
orbit's distances and phase angles, and moves them to r with band offsets
measured against Rubin r; objects without usable photometry get MPC H − 0.084
(the median difference). The details, including why the limit is 5.42, are
in `characterization/MPC/README.mpc`. `joint/coldfree_all_mpc.toml` adds the
MPC to the joint cold fit. `colour.toml` corrects V to Jester et al. (2005)
(V − r = +0.32 for cold colours).

## Cold classicals: free inclination

survey_debias package with the Survey Simulator container on CANFAR and provides routines to 'debias' the detected sample of survey and to code to find the MLE/MCMC constrained parameters of a Tapered Exponential LF that match the data from one or more surveys.  The joint fits in `joint/` select cold classicals one of two ways:

| Configurations | eps tag | Cold sample |
|---|---|---|
| `cold_*.toml` | `cold` | main-belt classical (`cl = cla`, `p = m`), ecliptic i < 5° |
| `coldfree_all.toml`, `coldfree_all_mpc.toml` | `coldfree` | main-belt classical, free inclination i_free < 4° |

The free inclination is that of Huang, Gladman & Volk (2022, ApJS 259, 54):
the angle to the local forcing pole, which a non-resonant orbit keeps while
its ecliptic i cycles, so the 4° cut separates the cold and hot components
much more cleanly. `survey_debias.ifree` gives it: Huang et al.'s Table 1
where the object is listed (540 of the 549 OSSOS++ main-belt classicals,
matched by name), otherwise their forcing-pole grid on barycentric elements.
JWST Sample A and N26 keep their papers' cold classifications (short arcs,
no usable node). The `coldfree` eps tables are planted with each model
row's free inclination kept (`survey_debias_eff --ifree-max 4`), because the
OSSOS 1.0 cold model is drawn about a forced plane, not the ecliptic.

## DEEP B1 (linked sample)

| | |
|---|---|
| Field | 10 DECam fields (B1a–j) along the invariable plane near RA 353°, Dec −4° |
| Epochs | 29 long stares, 2019-08-27 to 2021-10-06 (DEEP VI Table 1; times from the NOIRLab archive) |
| Detection | found on at least 4 stares (`detections_required = 4`) |
| Footprint | DECam CCD-row outline, 3.130 deg², fill 0.860 (60 searched CCDs) |
| Rates | 1.473–5.126 ″/hr, within 45° of the direction of decreasing ecliptic longitude |
| Efficiency per stare | DEEP III Eq. (3): c = 0.787, m25 = 26.22, k1 = 2.14, k2 = 5.02 (`function= double`) |
| Detections | 110: 73 cold, 23 hot, 6 resonant, 6 detached, 2 scattering |

VR is calibrated to Pan-STARRS r with no colour term, so the filter letter is
`r`. The efficiency is the survey-wide curve of DEEP III Figure 2 on every
night until DEEP's per-pointing-group table is available. DEEP needs an
`ossssim` with the equinox observer-frame fix and a `survey_debias` that aims
plants at every epoch of a K-of-N survey (both 2026-10-20). Sources, checks
and approximations are in `characterization/DEEP/README.deep`.

## Tests

```bash
python -m unittest discover -s tests
```

The simulator tests are skipped when the installed `ossssim` cannot load
multi-epoch survey roots.

## Licence

Released under the European Union Public Licence (EUPL), same as
SurveySimulator. See `LICENCE.txt`.
