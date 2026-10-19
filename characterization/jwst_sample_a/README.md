# JWST GO 1568 Sample A

Characterization of Sample A of Eduardo et al. (2026, AJ 172, 187; "E26"),
the 20 TNOs found in all three epochs of the GO 1568 NIRCam F150W2 search.
Morgan et al. (2026, AJ 172, 188) is the HST colour follow-up of the same
objects; it adds no detection samples or efficiencies, and it agrees with E26
on which objects are cold classicals.

```text
jwst_sample_a/
  README.md             this file
  jwst_sample_a.toml    survey_debias configuration
  jwst_sample_a.csv     detections: E26 Table 1 sample A
  survey.conf           detections_required = 3
  footprints.csv        MAST footprints of the 960 F150W2 SW exposures
  make_pointings.py     footprints.csv -> epoch{1,2,3}/pointings.list
  epoch1/ epoch2/ epoch3/
    pointings.list      one polygon: the three NIRCam strips, fill 0.91
    jwst_sample_a.eff   single-epoch efficiency (erfc lookup), filter= H
    JWST.csv            JPL Horizons barycentric vectors of JWST
```

## Where the values come from

| | Value | Source |
|---|---|---|
| Detections | 20 objects, F150W2 (AB) mean magnitudes | E26 Table 1, sample flag "A" (Section 5.5) |
| Epochs | JD 2459969.32121, 2459973.96788, 2459979.90857 | middle of each 20-tile visit sequence (MAST exposure times) |
| Footprint | 178.9 arcmin², three strips of 2.19′, 4.16′, 2.19′ by 20.9′ at PA 157.2° | MAST level-2 `s_region` of every SW detector |
| Fill factor | 0.91 | fraction of the strips imaged by one detector in both dithers, >15 px from its edge |
| Single-epoch efficiency | p0/2 erfc((m − 28.92)/0.61), p0 = 0.96 | E26 Eq. (8), Section 4.3.3 |
| Detection | all three epochs, so efficiency³ | E26 Section 5.5 and Figure 10 |
| Rate cut | 0.03–2.8″/hr (epoch 1), 0.03–1.4″/hr (epochs 2, 3) | KBMOD grid, E26 Section 4.1 |
| Photometric error | 0.11 × 10^(0.4(m − 28.6)) per epoch, catalogue = 3-epoch mean | fit to E26 Table 1 per-epoch errors |
| Magnitude limit | 30.5 | faintest implant, E26 Section 3.2 |

F150W2 has the filter letter `H`; `colour.toml` sets r − H = 1.0 (E26 Section
7.1, r − F150W2 ≈ 1.0).

## What was checked and what changed

### The sample is right; two classes were wrong

`jwst_sample_a.csv` lists exactly the objects flagged "A" in Table 1:
the 18 three-epoch TNOs found through the normal procedures (Sections 4 and
5.1) other than the targeted 2015 GK56 (JPB23), plus JPB01 and JPB05, which
the ResNet rejected in one epoch at S/N > 30 (Appendix B). JPB02 and JPB03
(linked only by bypassing the ResNet), JPB23, and the seven two-epoch objects
of sample B are correctly absent. Magnitudes, d_bary, i, a_min and e_min agree
with Table 1. The a and e columns are the minimum values for zero
line-of-sight velocity (Equations 14 and 15), not fitted elements.

Two `comp` values were changed:

- **JPB08** was `hot`; it is `cold`. It has d = 50.69 au and i = 3.33°,
  inside the E26 cold definition (38 < d < 55 au, i ≤ 5°, Section 7.1), and
  Table 1, Table 7 and Morgan et al. Table 1 all mark it cold.
- **JPB09** was `hot`; it is `centaur`. At d = 28.26 au it is the centaur
  that E26 leaves out of every luminosity-function sample (d > 30 au,
  Section 7.1; Table 7 "CC? No"). It stays in the file because it is a
  Sample A detection; select `comp` in {cold, hot} to reproduce the E26 fits.
  `centaur` uses the default colour group.

This gives 12 cold, 7 hot and 1 centaur.

There is no `.detections` file for this survey: the detection list is
`jwst_sample_a.csv`, and `survey_debias` writes
`jwst_sample_a.detections-full` from it.

### The efficiency parameters are Sample A's, but the functional form was wrong

The values 0.96, 28.92 and 0.61 are the single-epoch probability that an
implant on a detector in both dithers yields a good track (Eq. 8). They are
not a combined A/B efficiency. Sample A needs a good track in each of the
three epochs, and E26 show that the three-epoch recovery of implants is the
cube of Eq. (8) (Figure 10, open squares). Linking and the orbit-fit cuts
(min S/N > 5, Δm < 0.4, χ²_r < 25) kept all 413 correctly linked implant
triplets and lost 2 of 417 in linking, so they add nothing measurable. With
`detections_required = 3` the simulator draws the efficiency independently at
each epoch, so it applies exactly this cube. The 40% point of p_det³ is
m = 28.6, E26's quoted Sample A limit (Section 10). The 28.8 in the abstract
is the limit for samples A+B, and 28.92 is the single-epoch 48% point.

Eq. (8) is an **erfc**, not a tanh. The old `.eff` used `function= single`
(OSSSSim's (A/2)(1 − tanh((m − m0)/w))) with the same numbers. That makes the
transition too wide: the three-epoch efficiency was 9% too low at m = 28.0,
13% too low at 28.3 and 11% too low at 28.6, and too high faintward of 29.0
(×1.6 at 29.3). OSSSSim has no erfc, so the `.eff` now tabulates Eq. (8) in a
29-point lookup table. Interpolated and cubed, the table matches the formula to
better than 0.003.

The faintest Sample A object (JPB07, 28.60) is where you would expect it to
be. For the E26 slope α = 0.27, p_det³ predicts that 8% of detections, about
1.6 of 20, are fainter than 28.60, so seeing none has a probability of about
0.2. The single-epoch curve without cubing would predict 5 (P = 0.007), which
is where the impression of a sharp cutoff comes from.

### Same efficiency at every epoch; the KBMOD rate grid differs

E26 measure one Eq. (8) from implants in all epochs and find the three-epoch
recovery consistent with three independent, identical epochs. Lightcurve
amplitude has no measurable effect (Section 4.4). The three `.eff` files are
therefore identical except for `rate_cut`. The shift-and-stack grid was
−2.6 ≤ v_α ≤ 2.6, −1.3 ≤ v_β ≤ 1.3″/hr in epoch 1 and ±1.3″/hr in both
components in epochs 2 and 3 (Section 4.1). Along the strips, the direction
TNOs move, those boxes end at 2.8 and 1.4″/hr. The old cut, 0.03–8.66″/hr,
was the range of the implants, not the range searched. Only bodies inside
about 15 au move faster than 1.4″/hr, so the change does not affect any Sample A
object (all move at about 0.2 to 0.7″/hr).

The search thresholds do differ between the E26 samples:

- **Sample A:** S/N > 2.5 in each dither and χ²(f) < 7 per epoch, ResNet
  p > 0.65, then min S/N > 5 over the three-epoch orbit fit. The efficiency is
  A3 p_det³ with p_det from Eq. (8).
- **Sample B:** two epochs at S/N > 8. This has its own, brighter and sharper
  efficiency, p8 = 0.96, 28.74, 0.40 (Section 4.3.3), and two geometric
  channels, A2 p8² + A3 p8² (1 − p_det) (Section 5.5). A `jwst_sample_b`
  characterization needs those and cannot reuse these `.eff` files.
- **There is no sample C.** The "C" in the Table 1 samples column marks a
  dynamically cold object.

### Footprint: the real mosaic, not a 0.05 deg² square

The old `pointings.list` was a 13.4′ × 13.4′ square on the field centre with
fill 1. The mosaic is three strips about 21′ long, parallel to the ecliptic
(E26 Figure 1), separated by the two 40″ module gaps. About a third of it lay
outside the square, and the square included sky that was never imaged.
Because a moving object must be on a detector at all three epochs, the shape
matters: objects that move along a gap stay in it, and hot objects, which
drift across the strips, are lost more often than cold ones.

`make_pointings.py` builds the pointings from the MAST footprints
(`Mast.Jwst.Filtered.Nircam`, program 1568, F150W2, product level 2b, the
960 `*_cal.fits` rows in `footprints.csv`):

1. For each epoch, tile and SW detector, it takes the sky imaged in both
   dithers and removes a 15-pixel edge band (E26 count only apparitions more
   than 15 pixels inside a detector).
2. It fits the three strips to those 160 polygons and writes their outline as
   one S-shaped polygon, joined across the module gaps by 2″ bridges.
3. The fill factor is the fraction of the outline covered by the detector
   polygons: 0.905, written as 0.91. That fraction is lost to the 4″ chip
   gaps, the edge bands and the dither offsets.

One polygon per epoch is deliberate. `survey_debias` aims plants uniformly
over each tile's bounding box and weights them by the summed box area, which
is only unbiased when tiles do not overlap. The 160 detector polygons per epoch
overlap at tile edges.

The three-epoch geometric area was checked by planting circular orbits at
35–55 au uniformly over the field and running them through `ossssim` with
unit efficiency:

| Three-epoch area A3 [arcmin²] | Cold (i < 5°) | Hot (i > 5°) |
|---|---|---|
| E26 Section 5.5, from implants | 139 | 128 |
| Old square, fill 1 | 164 | 163 |
| 160 detector polygons per epoch, each at its visit time | 115 ± 2 | 95 ± 2 |
| S-polygon, fill 1 | 148 | 128 |
| S-polygon, fill 0.91 (this characterization) | 115 ± 2 | 98 ± 2 |

With fill 0.91 the S-polygon reproduces the detector-level geometry to 0.3%
(cold) and 3% (hot). The best three-epoch fill is 0.920 (cold) or 0.906
(hot). The old square overstated the Sample A area by 42% (cold) and 72% (hot).
Fitting a fill factor to the old square cannot fix this, because the square
gives the same area for cold and hot objects. Computing a single time for each
epoch, rather than one per tile, changes A3 by 2%.

In the same simulation, the ratios A1 : A2 : A3 agree with E26's (cold
0.12 : 0.28 : 0.59 against 0.13 : 0.25 : 0.62). The absolute areas, however,
are about 20% (cold) and 25% (hot) below E26's. E26's areas are what this
footprint gives when the chip gaps count as imaged (the S-polygon with fill 1:
148 and 128). For none of the motions tried does the real footprint give both
A3 = 139 and A1 + A2 + A3 = 224. This characterization follows the footprint, because E26's p_det is
conditional on the source landing on a valid part of a detector. To reproduce
E26's normalization instead, set the fill factor to 1.0, which gives about 6%
more than 139 for cold objects and matches 128 for hot ones. Only the
luminosity-function normalization depends on this choice; the slope does not.

### Photometric errors and the catalogue magnitude

The old `mag_error` (0.05 at 28.2 rising to 0.5 at 29.5) is not in E26. The
per-epoch uncertainties in Table 1 follow the background-limited form
σ(m) = 0.11 × 10^(0.4(m − 28.6)), with 0.035 dex rms: 0.02 at 26.7, 0.04 at
27.5, 0.11 at 28.6 and 0.21 at 29.3. The old curve was a factor of 3 too small
at 27.5 and a factor of 2 too large at 29.5. The Table 1 magnitude is the mean
over all three epochs, so `phot_frac= 0.0 0.0 1.0` makes the simulator report
the mean of three draws (σ/√3, e.g. 0.06 at 28.6, as for JPB06 and JPB07).

### End-to-end check

On the line of sight through the field centre (44 au), the fraction of
simulated objects found is 0.648, 0.629, 0.545, 0.317, 0.098 and 0.013 at
m = 25.0, 28.0, 28.3, 28.6, 28.9 and 29.2, with 1000 draws each. The expected
values, (0.91 p_det(m))³, are 0.667, 0.634, 0.527, 0.305, 0.093 and 0.011.
`survey_debias` loads the three epochs with selection S(m) = p_det³, for
example 0.405 at 28.6.

## Not modelled

- The slight drop in Figure 10 at the bright end of the implant range, which
  E26 attribute to the PSF model. The efficiency is flat at 0.96 brightward of
  27.3. The brightest Sample A object is 24.47.
- The extra chance of detection for an object imaged on two overlapping tiles
  in one epoch (E26 found four such TNOs). The S-polygon counts the overlaps
  once.
- JPB01 and JPB05 are in the sample although the ResNet rejected one of their
  epochs. E26 argue that the implant efficiency still applies to them.
