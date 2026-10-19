# MA: Alexandersen et al. (2016) CFHT survey

Characterization of the low-lat (`mal`) and high-lat (`mah`) blocks of
Alexandersen et al. (2016, AJ 152, 111; "A16"; arXiv:1411.7953), a 32 deg²
CFHT/MegaCam r-band survey near Neptune's L4 point, designed for plutinos and
Neptunian trojans.

```text
MA/
  README.md                       this file
  pointings.list                  4 low-lat + 2 high-lat rectangles at the discovery epochs
  ma11la.eff                      efficiency, rate cut and tracking of the low-lat block
  ma12ha.eff                      efficiency, rate cut and tracking of the high-lat block
  MA.detections                   the 77 characterized objects (release format)
  Alexandersen_2016_AJ_152_111.pdf
  ajaa34e0t3_mrt.txt              A16 Table 3, machine-readable: all 87 detections
  ajaa34e0survey/                 A16 supplementary survey package
    MA11l12h/                     the characterization used in A16
      pointings.list, ma11la.eff, ma12ha.eff, README
    Res31Model.f, Res32Model.f, Res41Model.f, CoorbModel.f
                                  model generators for the 3:1, 3:2, 4:1 and co-orbital runs
    param*.txt                    their parameter files
    KRQ11hot_resaei_100Myr_incfix.out
                                  scattering input to CoorbModel.f (Alexandersen et al. 2013)
    runme.sh                      example survey-simulator run
    ReadMe
  ajaa34e0astrometric/            A16 supplementary astrometric catalogue
    LL_cur.*, HH_cur.*            low-lat and high-lat block catalogues and exposure lists
    ReadMe                        (describes the survey package, not the catalogue)
```

## Sources

The supplementary files are the source of record. `ajaa34e0survey/MA11l12h/`
holds the pointings and efficiency files used for the A16 simulations, and
`ajaa34e0t3_mrt.txt` is Table 3 at full precision with the
characterization flag, rate, phase angle and classification of every object.
The printed paper (Tables 1–3, Sections 2–3) agrees with them, except that
the printed Table 2 efficiency fits are not the ones in the supplementary
`.eff` files (see below). The ReadMe of the survey package lists
`Res11Model.f` and `param11.txt`, but they are not in the package.

See `CFEPS_g/README.formats` for the block format.

## The survey (A16 Sections 2 and 3)

| | Low-lat block | High-lat block |
|---|---|---|
| Area | 20 deg² (5° × 4°) near the ecliptic | 12 deg² (4° × 3°), about 15° north |
| Sub-blocks | south half, north half (10 deg² each) | west half, east half (6 deg² each) |
| Discovery nights | 2011 Oct 24 (south), Oct 26 (north) | 2012 Oct 20 (west), Oct 21 (east) |
| Image quality | mostly 0.5–0.6″ | 0.5–0.8″ |
| Detected / characterized | 70 / 62 | 17 / 15 |

- **Discovery images:** triplets of 320 s r-band exposures about 1 hr apart.
  The slowest detectable motion is about 0.3″/hr. High-lat triplets are
  interleaved with nailing images of the other sub-block.
- **Pipeline:** the CFEPS moving-object pipeline (Petit et al. 2004). Both
  detection algorithms must find a candidate, then every candidate is vetted
  by eye.
- **Implants:** artificial objects went into time-scrambled copies of the
  triplets, not into the searched images. Real and fake candidates were vetted
  blind and mixed together. Of 21,858 accepted fake candidates, 3 were false
  positives, all fainter than the characterization limit, so the
  characterized sample has zero false positives.
- **Characterization limit:** for each rate range, the magnitude at which the
  efficiency is 40% of its maximum, rounded up to 0.01 mag.
- **Tracking:** in the discovery year each whole block was shifted across the
  sky at the mean TNO rate, so tracking does not depend on assumed orbits. The
  high-lat block was also imaged in 2011 Aug, 2011 Nov and 2012 Jan, which
  serve as precovery. Arcs reach 28 months. Of the 77 characterized objects,
  two were not tracked:
  - **mal11nt** (2011 UU412, m_r = 24.5) is near the limit and was invisible
    in the shallower tracking images.
  - **mah11nt** (m_r = 24.2) sheared off the field during the discovery year.
- **Photometry:** m_r is the mean r magnitude from the moving-object pipeline
  over the discovery triplet, and H_r includes the phase correction. Where A16
  combine with CFEPS they assume g − r = 0.5 for plutinos.

## What was checked and what changed

### Efficiencies: the supplementary files, with a key typo fixed

`ma11la.eff` and `ma12ha.eff` are now the supplementary
`MA11l12h/ma11la.eff` and `ma12ha.eff`, with one change. In `ma11la.eff`
the last two rate ranges were written `rates:` instead of `rates=`. The
simulator (the F77 reader of the time and the current F95 one alike) skips
any line without `=`. Taken literally, therefore:

- the 5.38–6.40 and then the 6.40–10.00 fits overwrote the 3.72–5.38 range,
  which ended up with f21 = 0.786, m_e = 24.375 and limit 24.37 instead of
  0.901, 24.560 and 24.57;
- rates from 5.38 to 10″/hr had no efficiency function, so efficiency 0.
  The Uranian trojan mal01 (6.27″/hr) could not have been simulated as
  detected.

The two lines are written `rates=` here, so all five ranges load. The
supplementary `ma12ha.eff` is used unchanged.

| Key | Supplement (now used) | Previous file |
|---|---|---|
| `phot_frac` (both) | `1.000 0.000 0.000` | `0.026 0.156 0.818` |
| `track_frac` low | `1.000 24.52 -3.333` | `1.000 24.500 -2.78` |
| `track_frac` high | `0.933 24.70 -5.000` | `0.930 24.600 0.0` |
| rate ranges low | 0.50–2.06, 2.06–3.72, 3.72–5.38, 5.38–6.40, 6.40–10.00 | Table 2: …, 5.38–10.00 |
| rate ranges high | 0.50–1.00, 1.00–6.00, 6.00–10.36 | Table 2: 0.50–2.06, 2.06–6.00, 6.00–10.36 |

`rate_cut`, `mag_error` and `filter` are the same in both. `track_frac`
gives min(track_max, 1 + (m − track_mag) × slope), so the supplement's
tracking is 100% to m_r = 24.52 and then falls to 0.93 at 24.54 and 0.67 at
the 24.62 limit (low-lat), and a constant 93.3% up to 24.71 (high-lat). This
is A16 Section 3.3. A16 say the triplet mean defines m_r, while
`phot_frac= 1 0 0` simulates a single measurement; the supplement's value is
kept.

The printed Table 2 rounds the same fits in the ranges that coincide. In the
ranges that do not coincide, it gives fits that are not in any supplementary
file:

| Block | Rates [″/hr] | `.eff` (supplement): f21, k, m_e, w, limit | A16 Table 2 |
|---|---|---|---|
| low | 0.50–2.06 | 0.950, 0.0136, 24.620, 0.126, 24.62 | 0.96, 0.0148, 24.62, 0.125, 24.62 |
| low | 2.06–3.72 | 0.921, 0.0126, 24.611, 0.137, 24.61 | 0.91, 0.0118, 24.60, 0.133, 24.61 |
| low | 3.72–5.38 | 0.901, 0.0103, 24.560, 0.147, 24.57 | 0.90, 0.0106, 24.55, 0.142, 24.57 |
| low | 5.38–6.40 | 0.994, 0.0481, 24.599, 0.067, 24.42 | 5.38–10.00: 0.88, 0.0136, 24.48, 0.130, 24.49 |
| low | 6.40–10.00 | 0.786, 0.0135, 24.375, 0.137, 24.37 | (in 5.38–10.00) |
| high | 0.50–1.00 | 1.002, 0.0355, 24.677, 0.080, 24.59 | 0.50–2.06: 0.92, 0.0194, 24.64, 0.118, 24.60 |
| high | 1.00–6.00 | 0.911, 0.0146, 24.608, 0.120, 24.60 | 2.06–6.00: 0.90, 0.0133, 24.59, 0.115, 24.60 |
| high | 6.00–10.36 | 0.881, 0.0136, 24.482, 0.130, 24.48 | 0.88, 0.0139, 24.49, 0.124, 24.49 |

Where they differ, the supplement's efficiency is 10–30% lower than Table 2
at m = 23.5–24.3 for low-lat rates above 5.38″/hr. That is, objects inside
about 25 au, such as mal01 at 20.3 au and 6.3″/hr. For high-lat rates below
1″/hr (objects beyond about 150 au) it is 8–9% lower at m = 24.0–24.3. For
1–2″/hr it is 3–4% higher at m = 24.0–24.3.

### Pointings: low-lat from the supplement, high-lat recentred

**Low-lat.** The supplement (and Table 1) lists the low-lat block as four
0.995°-high rows. `pointings.list` now uses them as written. Before, the
rows were merged in pairs into two 1.990°-high rectangles, which
gave the same area to 0.1%. The CADC MegaPrime archive confirms the block: 20 fields at
RA 25.8–29.8° and Dec +13.4° to +16.4°, the southern ten with triplets on
2011 Oct 24 and the northern ten on Oct 26 (six of them also on Oct 25).

**High-lat.** The supplement (and Table 1) puts the two halves at 1.993° ×
2.984°, centred at 01:41:01.70 and 01:49:58.30, +28:10:00. Five of the 15
`mah` objects (mah01, mah02, mah05, mah10, mah12) were discovered at
RA 23.2–23.8°, west of both rectangles. Horizons confirms the catalogue
positions (2012 UV177 = mah02 at RA 23.28°, Dec +28.68° on the discovery
night). The CADC archive has the discovery block as 12 r-band 320 s fields
on a 4 × 3 grid at RA 23.6, 24.7, 25.8, 26.9° and Dec +27.0, +28.0, +29.0°.
The western six have their triplet on MJD 56220 (mean 56220.404) and the
eastern six on MJD 56221 (mean 56221.451), which matches the supplement's
JDs. The extent of each half, from the CCD polygons of its first triplet
exposures:

| Half | CCD footprint | Supplement / Table 1 |
|---|---|---|
| West | 1.983° × 2.996° at 01:36:35.9 +28:01:04 | 1.993° × 2.984° at 01:41:01.7 +28:10:00 |
| East | 1.983° × 2.996° at 01:45:26.8 +28:01:04 | 1.993° × 2.984° at 01:49:58.3 +28:10:00 |

The same measurement reproduces the low-lat rows to 0.01° in size and 20″ in
position. The supplement's high-lat centres are therefore one MegaCam field
(0.98°) too far east and 9′ too far north. `pointings.list` uses the
measured centres with the supplement's sizes, JDs and fill factor (0.907). The
CCD polygons give a fill of about 0.95. The supplement's lower value
presumably also excludes masked and unsearched pixels, so it is kept. All 15
`mah` objects are now inside the block. In
`survey_debias_eff --check-detections` they are recovered at the expected
rate (efficiency × 0.933 tracking × 0.907 fill). With the supplement's
centres, five could never be recovered.

### Detections: release values kept, mal58 reclassified

`MA.detections` is in the OSSOS release format. It holds exactly the 77
objects with characterization flag 1 in `ajaa34e0t3_mrt.txt` (62 `mal`, 15
`mah`, including mal11nt and mah11nt), and none of the ten with flag 0.
Against the machine-readable Table 3:

- **Apparent magnitudes** and **discovery distances** agree. **Rates**
  agree to 0.01–0.02″/hr.
- **e_mag** is the supplement's e_rmag divided by √2 in every row. The
  simulator does not read it (photometric errors come from `mag_error`).
  Left as is.
- **Hsur** is 0.01–0.04 mag brighter than the supplement's phase-corrected
  rMag. That follows the release's convention for Hsur, so it is left as is.
- **Orbits:** five objects have release arcs of about 4 yr, against 1.3–2.2 yr
  in the supplement: mal17, mal23, mal35, mal47 and mal53. Their elements
  differ accordingly, most for mal23 (a = 61.84 against 62.12 au, e = 0.433
  against 0.438). Eleven others differ in a by at most 0.0003 au. The
  release orbits use the longer arcs and are kept.
- **Classes:** the release classifications differ from A16 for eight
  objects. Only one changes the cold-classical sample (main classical,
  i < 5°), so only that one was changed:

| Object | a, e, i | Release (`cl p j k sh`) | A16 | Changed? |
|---|---|---|---|---|
| mal58 | 47.88, 0.161, 2.61° | `cla m` insecure | outer classical, secure | yes, to `cla o … S` |
| mal01 | 19.09, 0.177, 10.81° | `cen` | 1:1 Uranus, secure | no |
| mah01 | 30.07, 0.259, 53.89° | `cen` | 1:1 Neptune (temporary), secure | no |
| mah05 | 40.67, 0.117, 19.24° | 11:7 insecure | main classical, secure | no |
| mah09 | 42.82, 0.099, 16.60° | 17:10 insecure | main classical, secure | no |
| mah08 | 92.32, 0.617, 24.85° | detached | 16:3 insecure | no |
| mah14 | 73.79, 0.492, 16.35° | 23:6 insecure | detached, insecure | no |
| mal19 | 67.32, 0.483, 16.71° | 10:3 insecure | detached, secure | no |

mal58 lies outside the 2:1 resonance (47.7 au), where the Gladman et al.
(2008) scheme that A16 use, and the OSSOS release itself, call a classical
"outer". The release `m` is inconsistent with its own a, and it put mal58
in the `cla`/`m`/i < 5° cold sample used in `joint/cold_all.toml`. With the
change, MA has ten cold classicals (mal36, 39, 42, 43, 46, 48, 51, 52, 55, 56),
as in A16. The other seven differences are naming conventions (the release
calls a < a_Neptune objects centaurs) or reclassifications of insecure
resonances. All have i > 10°, so they do not affect the cold sample and are
left as in the release.

## Not modelled or still to do

- The direction cone of `rate_cut` (−21.05 ± 15.62° and −20.655 ± 19.335°
  from due west) is not described in A16. The detected objects move at
  −19° (low-lat) and −21° (high-lat), well inside both cones.
- `eps_cold*.csv` and `plants_cold.npz` were made with the earlier pointings
  and efficiencies and need to be regenerated, and the joint fit rerun (mal58
  has left the cold sample).
