# Crystalline silicon

Crystalline silicon (`Si`) has one recommended model of its complex refractive index,
from the extreme ultraviolet to the far infrared, at 20 degrees C. This page gives:

- the model and what it rests on;
- how well it is known in each range;
- how it compares with the published alternatives.

The manual for solids is {doc}`solid`. The register entries are TC-1075 to TC-1087
(the sources) and TC-1088 (the fit).

```python
from optmat.solid.models.registry import default_model, models

silicon = default_model("Si")  # fit/Si, at 20 degrees C
n = silicon.refractive_index(1550e-9)  # complex n + ik
beta2 = silicon.gvd(2000e-9)  # s^2/m

cold = models("Si")["frey_2006/sellmeier_Si"].at_temperature(120.0)  # 20-300 K
```

## The complex index

**The recommended model** is `fit/Si` (TC-1088). It is Franta et al.'s universal
dispersion model of float-zone silicon, refitted to the measurements, in the same way
as {doc}`silica`'s `fit/SiO2`.
- **Range:** 31 nm to 143 µm, at 20 degrees C.
- **Physics:** it is causal, passive and normalised to its transition strengths, with
  an exact d eps/d omega.
- **The start:** the newAD2 parameter file of float-zone c-Si
  (`franta_2017/udm_c_si`), which the authors' site describes as unpublished. optmat
  reproduces the authors' own table of it to the last printed digit (TC-1075). It is
  0.016 above the absolute refractometry in the near infrared, and its k at 633 nm is
  2.6 times the measured.
- **What it rests on:**
  - in the transparent window (1.2-5.5 µm), the absolute indices of Frey, Leviton &
    Madison (CHARMS, ±1e-4), taken from 295 K to 293.15 K with their own measured
    dn/dT;
  - on the absorption edge (250-1200 nm), Schinke et al.'s absorption coefficient,
    whose uncertainty analysis gives 0.2-13 % (k = 2);
  - elsewhere, Franta et al.'s temperature-dependent model at 20 degrees C, the
    numbers the authors gave to the refractiveindex.info database (TC-1075). Its
    ultraviolet comes from their ellipsometry and its vacuum ultraviolet (to 40 eV)
    from their synchrotron reflectance. Its infrared k is their measured lattice
    absorption. Its window n is tied to Frey et al.'s and is not fitted beside it.
- **What it changed:** the interband, higher-energy and Urbach strengths, the
  localised-state strengths at the E1 and E2 critical points, the four excitons, and
  one factor on the phonon bands (x1.006).

<!-- overview:start -->
| lambda (nm) | n | k |
|---|---|---|
| 50 | 0.7456 | 0.0542 |
| 100 | 0.2444 | 0.970 |
| 150 | 0.5827 | 1.81 |
| 200 | 0.9932 | 2.76 |
| 250 | 1.6372 | 3.62 |
| 300 | 5.0216 | 4.12 |
| 370 | 6.8785 | 2.10 |
| 400 | 5.6010 | 0.331 |
| 500 | 4.2961 | 0.0474 |
| 633 | 3.8704 | 0.0180 |
| 800 | 3.6776 | 0.00598 |
| 1000 | 3.5744 | 0.000491 |
| 1064 | 3.5540 | 8.99e-05 |
| 1550 | 3.4783 | 6.43e-10 |
| 2000 | 3.4532 | 1.40e-10 |
| 5000 | 3.4226 | 1.72e-07 |
| 10000 | 3.4182 | 6.05e-05 |
| 20000 | 3.4176 | 0.000268 |
| 100000 | 3.4173 | 0.000202 |
<!-- overview:end -->

## The transparent window

The absolute refractometry and the models at 20 degrees C:

<!-- window:start -->
| | Frey 2006, table | Frey 2006, formula | Li 1980 | Franta 20 C | newAD2 | fit/Si |
|---|---|---|---|---|---|---|
| n at 1.5 µm | 3.48268 | 3.48281 | 3.47990 | 3.48262 | 3.49947 | 3.48269 |
| n at 2.0 µm | 3.45321 | 3.45320 | 3.45100 | 3.45323 | 3.46896 | 3.45323 |
| GVD at 2.0 µm (fs^2/mm) | | 817.9 | | | 844.7 | 814.8 |
<!-- window:end -->

- **Frey et al.** measured in vacuum, from 30 to 295 K and 1.1 to 5.5 µm (TC-1079).
  The columns are their 295 K table taken to 20 degrees C with their dn/dT (Table 4,
  -3.3e-4 at 2 µm), and their temperature-dependent Sellmeier formula
  (`frey_2006/sellmeier_Si`). The formula reproduces their own tables only to about
  1e-4 (rms 1.35e-4, up to 1e-3 at 1.1 µm), so the fit uses the tables.
- **The fit** follows Frey's tables to within 1.4e-4 from 1.2 to 5.5 µm. The largest
  departure is at 1.2 µm, on the absorption edge, and it is +0.6-1.2e-4 at 2.5-5.5 µm,
  where the model's shape does not quite follow the data.
- **Li's** critical review (1980) lies 2-4e-3 below Frey's, beyond its own stated
  ±2e-3 (TC-1080). It is a check, not fitted. Li does not say whether his index is
  absolute or relative to air; read as relative to air, it is still 1.3-3.3e-3 low.
- **Green's** self-consistent table (300 K, 2008) lies 1.7-3.3e-3 below the fit at
  0.93-1.45 µm. From Frey's dn/dT (1.8e-4 per K) it should lie 1.2e-3 above at 300 K,
  so it is 3-4.5e-3 below Frey's level (TC-1081).
- **GVD:** the fit's GVD is 1112 fs^2/mm at 1550 nm and 815 at 2 µm. That is within
  4 % of Frey's formula from 1.55 to 5 µm, and 6 % above it at 1.2 µm. The newAD2
  file is 3-6 % higher.
- **Temperature:** dn/dT is about 1.8e-4 per K at 2 µm, so 1 K moves n by twice
  Frey's uncertainty. For other temperatures in the window, use Frey's formula
  `frey_2006/sellmeier_Si` with `at_temperature` (20-300 K, 1.1-5.6 µm). It carries
  the formula's own 1e-4 misfit. Their dn/dT at all twelve temperatures is curated as
  datasets.

## The absorption edge

Silicon's gap is indirect, at 1.12 eV (1.1 µm). Below the direct transitions at
3.4 eV (365 nm), its absorption rises by six decades through phonon-assisted
transitions.

<!-- edge:start -->
| lambda (nm) | Schinke 2015 (295 K) | Green 2008 (300 K) | Franta 20 C | newAD2 | fit/Si |
|---|---|---|---|---|---|
| 400 | 0.326 | 0.296 | 0.327 | 0.335 | 0.331 |
| 500 | 0.0485 | 0.0442 | 0.0523 | 0.0693 | 0.0474 |
| 633 | 0.0158 | 0.0161 | 0.0152 | 0.0398 | 0.0180 |
| 800 | 0.00527 | 0.00541 | 0.00598 | 0.0166 | 0.00598 |
| 900 | 0.00211 | 0.00217 | 0.00260 | 0.00702 | 0.00218 |
| 1000 | 0.000490 | 0.000509 | 0.000464 | 0.00143 | 0.000491 |
| 1100 | 3.02e-05 | 3.06e-05 | 3.05e-05 | 3.44e-05 | 3.04e-05 |
| 1200 | 1.39e-07 | 2.10e-07 | 1.16e-07 | 5.18e-07 | 3.48e-07 |
<!-- edge:end -->

The columns give k; the tables of alpha are converted with k = alpha lambda / (4 pi).

- **The fit's k** is within 4 % of Schinke's below 370 nm, within -8 to +13 % at
  370-550 nm and 810-1100 nm, 5-19 % high at 550-810 nm, and up to 27 % low at
  1100-1150 nm. At 550-810 nm Schinke's own uncertainty reaches 13 % (k = 2), and the
  values come from a Tauc fit. The model's few smooth terms cannot follow the shape
  more closely.
- **Green's** table is a Kramers-Kronig re-analysis of others' measurements,
  including Jellison's and Aspnes's (TC-1081). It is a check, so that those
  measurements are not counted twice. It agrees with Schinke's within 10 % to
  1.1 µm.
- **Beyond 1.15 µm** the band-to-band absorption falls by five decades in 0.15 eV:
  alpha < 0.1 cm^-1, k < 1e-7. The model's Urbach tail cannot follow it. It is 2.5
  times Schinke's at 1.2 µm and levels off at 1e-10 to 1e-9 in the window. **The
  model's k beyond 1.15 µm is an upper bound, not a measurement.** In a real wafer
  the free carriers of the doping usually absorb more.
- **Temperature:** k on the edge changes by 0.5-2 % per K at 0.9-1.15 µm (Green's
  normalised coefficients, 54-210e-4 per K). Schinke's 295 K data are fitted at 20 degrees C
  with a 1 % allowance.

## The ultraviolet

The ellipsometric dielectric functions bracket the model across the E1 (3.4 eV) and E2
(4.3 eV) critical points:

<!-- uv:start -->
| E (eV) | Jellison n | Jellison k | Aspnes n | Aspnes k | Franta n | Franta k | fit n | fit k |
|---|---|---|---|---|---|---|---|---|
| 3.0 | 5.226 | 0.196 | 5.222 | 0.269 | 5.239 | 0.222 | 5.243 | 0.219 |
| 3.4 | 6.486 | 2.674 | 6.522 | 2.705 | 6.490 | 2.712 | 6.480 | 2.727 |
| 4.0 | 5.077 | 3.557 | 5.010 | 3.586 | 5.027 | 3.571 | 5.014 | 3.577 |
| 4.3 | 4.191 | 5.280 | 4.087 | 5.395 | 4.095 | 5.303 | 4.084 | 5.305 |
| 5.0 | 1.683 | 3.608 | 1.570 | 3.565 | 1.634 | 3.548 | 1.627 | 3.555 |
<!-- uv:end -->

- **Aspnes & Studna's** (1983) n is 0.10-0.13 below Jellison's at 4.3-5.3 eV (TC-1083,
  TC-1084). Herzinger et al. (1998) found that Aspnes's table behaves as if 0.18 nm of
  oxide remained on the surface (TC-1085). Aspnes's k below 2.8 eV is Dash & Newman's
  absorption, not his measurement, and is left out.
- **Jellison's** (1992) error limits are ±0.012-0.046 in n. Two misprinted n values
  (595 and 705 nm) are corrected as errata.
- **The fit** follows Franta's ultraviolet n to 3e-3 on average above 120 nm and
  7e-3 at 31-120 nm, and to 0.04 at the sharpest structure. It lies between Jellison's and Aspnes's n, nearer Aspnes's at
  4-5 eV. In this range the spread of the measurements, 0.01-0.1 in n and a few per
  cent in k, is the uncertainty.
- **In the vacuum and extreme ultraviolet** (31-200 nm) Franta's model is fitted to
  their own synchrotron reflectance (6-40 eV). Edwards's chapter in Palik's handbook
  compiles Philipp's (1972) Kramers-Kronig analysis and Hunter's films to 100 eV
  (TC-1086). Philipp's n is 23 % below Franta's at 150 nm, and 4 % above at 210 nm.

## The infrared

From 5.6 to 143 µm the model's n is Franta's to 1.5e-4. Its k is Franta's measured
lattice absorption, mostly two-phonon bands, from 6 to 64 µm: to 10 % rms, and up to
45 % on the sharpest bands. These are
float-zone values. Czochralski silicon adds the 9 µm (1107 cm^-1) band of interstitial
oxygen, and doped silicon adds free-carrier absorption, both sample-dependent, which
the model does not include.

## Limitations

- 20 degrees C only. The only temperature-dependent model is Frey's formula, in the
  window.
- Franta et al.'s full temperature-dependent model (2017) follows the indirect edge's
  phonon structure from 10 to 700 K. It cannot be implemented: its parameters, more
  than a thousand, are not printed in the 2017, 2018 or 2019 papers (TC-1075,
  TC-1077, TC-1078). The refractiveindex.info tables at 10-500 K are its only
  numerical form.
- The model is intrinsic, float-zone silicon. It has no oxygen, no doping and no free
  carriers.
