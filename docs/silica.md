# Fused silica

Fused silica (silica glass, `SiO2`) has linear and nonlinear models. This page covers:

- its complex refractive index from the extreme ultraviolet to the far infrared;
- its nonlinear index n2 with the Raman part;
- how well both are known, and how they compare with the published alternatives.

The manual for solids is {doc}`solid`; the register entries are TC-321 to TC-324 and
TC-399 (the linear index) and TC-393 (n2).

```python
from optmat.solid.models.registry import default_model, models
from optmat.solid.kerr import kerr_models

silica = default_model("SiO2")  # fit/SiO2, at 20 degrees C
n = silica.refractive_index(800e-9)  # complex n + ik
beta2 = silica.gvd(800e-9)  # s^2/m
deps = silica.dpermittivity_domega(120e-9)  # exact d eps / d omega

cold = models("SiO2")["leviton_2008/sellmeier"].at_temperature(120.0)  # 30-300 K

kerr = kerr_models()["kerr/SiO2"]
kerr.n2(1030e-9)  # long-pulse n2, m^2/W
kerr.n2_effective(1030e-9, 50e-15)  # what a 50 fs pulse measures
```

## The complex index

**The recommended model** is `fit/SiO2` (TC-399). It starts from Franta et al.'s
universal dispersion model of bulk Schott Lithosil Q2 (`franta_2016/udm_lithosil_q2`),
refitted to the absolute refractometry of the transparent window.
- **Range:** 24.8 nm to 125 µm, at 20 degrees C.
- **Physics:** it is causal, passive and normalised to its transition strengths, with
  an exact d eps/d omega.
- **What it rests on:**
  - in the window, the absolute indices of Leviton & Frey (Corning 7980, 2008) and
    Malitson's three prisms (1965), which agree with each other to 1e-5;
  - Franta's measured absorption at the band edge (146-225 nm) and in the infrared,
    and his infrared n;
  - in the vacuum ultraviolet, Franta's model, which above 8.5 eV was fitted to
    Philipp's tables in Palik's handbook.
- **What it changed:**
  - the interband strength (-2 %) and the higher-energy term (-6 %, and 16.1 eV
    against 17.0 eV);
  - the band edge (8.355 against 8.363 eV) and the Urbach width (0.188 against
    0.189 eV);
  - one factor each on the excitons (x1.59), the phonon bands (x1.0006 and x1.003)
    and the multiphonon bands (x1.12).
- **What it keeps** of Franta's model (rms):
  - the band-edge k to 1.7 %;
  - the vacuum-ultraviolet n and k to 0.008 and 1.7 %;
  - in the infrared, k to 6-9 % and n to 2e-3.

**In the vacuum ultraviolet** the band edge (E_g = 8.36 eV, 148 nm, with an Urbach tail)
and the excitonic structure near 10-12 eV dominate:

<!-- vuv:start -->
| lambda (nm) | n | k |
|---|---|---|
| 50 | 0.8080 | 0.300 |
| 80 | 1.1935 | 0.650 |
| 100 | 1.4305 | 0.767 |
| 110 | 1.8046 | 0.739 |
| 120 | 2.2399 | 0.969 |
| 130 | 2.0202 | 0.0684 |
| 150 | 1.7061 | 0.000880 |
| 160 | 1.6490 | 6.68e-05 |
| 180 | 1.5859 | 9.12e-07 |
| 200 | 1.5511 | 2.98e-08 |
<!-- vuv:end -->

Extrapolating a Sellmeier formula through a pole does not describe absorption.
For example, evaluating Malitson's printed formula at vacuum wavelengths with a
complex square root, outside its stated range, gives Re n = 0.000 at 110 nm and 2.94 at 120 nm,
where the recommended complex model gives 1.80 and 2.24. The discrepancy is not an
uncertainty estimate for either model.

**In the transparent window** the absolute refractometry agrees to 1e-5, once each
measurement is referred to vacuum:

<!-- window:start -->
| At 800 nm | Leviton 2008 (20 C) | Malitson 1965 | Arosa 2020 | Franta 2016 | fit/SiO2 |
|---|---|---|---|---|---|
| n | 1.45371 | 1.45371 | 1.45337 | 1.45355 | 1.45371 |
| group index | 1.46757 | 1.46756 | 1.46727 | 1.46726 | 1.46755 |
| GVD (fs^2/mm) | 36.31 | 36.23 | 35.86 | 35.76 | 36.20 |
<!-- window:end -->

- **Relative to air or to vacuum.** Refractometry by minimum deviation in room air
  gives the index relative to that air, 2.7e-4 below the absolute index at room
  conditions (TC-398).
  - **Leviton & Frey** measured in vacuum (the CHARMS refractometer), from 30 to
    310 K, to +-1e-5 (`leviton_2008`).
  - **Malitson's** 1965 indices are relative to air. The paper does not say so, but
    Leviton & Frey read them so, and converted to vacuum they agree with Leviton's to
    4.5e-6 rms from 0.4 to 2 µm. `malitson_1965/sellmeier` now returns the absolute
    index, through Edlen's air at 20 degrees C (TC-399).
  - **Rodney & Spindler's** 1954 indices, from the same apparatus, are relative to air
    too (`rodney_1954`).
  - **Arosa & de la Fuente's** formula (2020) is levelled by an Abbe refractometer,
    which also reads relative to air. Taken to vacuum it agrees with Leviton's to
    5e-5 from 0.4 to 2 µm; the table shows it as published.
- **Franta's published model** lies 3.4e-4 below Leviton's at 400 nm and 0.8e-4 below
  at 2 µm, and its GVD is 1.5 % low at 800 nm. It was the recommended model until this
  comparison (TC-324); the refit keeps everything else about it.
- **Uncertainty of the recommended model in the window (0.21-2.6 µm):**
  - n to about 2e-5: the fit's residuals against Malitson and Leviton are 1-2e-5, and
    Malitson's twelve specimens of three grades scatter by 1-3e-5 about it;
  - the group index to about 3e-5 (the refit is 2.3e-5 below Leviton's formula and
    0.6e-5 below Malitson's at 800 nm);
  - GVD to about 0.3 %: Leviton's formula gives 36.31 fs^2/mm and Malitson's 36.23 at
    800 nm, against the fit's 36.20.
- **Near the 2.73 µm OH band and beyond** grades differ: the fit, with Lithosil Q2's
  OH band and Franta's infrared, departs from Malitson's glasses by up to 5e-5 at
  3.3 µm.
- **The dispersion near 2 µm carries Lithosil Q2's OH bands.** The fit keeps Franta's
  bands of the plate unscaled (Nph36-52; TC-399), and their Kramers-Kronig partners shape
  the dispersion well before the bands themselves:
  - the 2.73 µm fundamental (Nph36-41);
  - a narrow group at 2.21-2.25 µm, 51-293 cm⁻¹ wide, with k up to 1.0e-5 (0.6 cm⁻¹) at
    2212 nm (Nph45-47);
  - weaker bands at 1.90-2.59 µm and 1.38 µm.

  At 2 µm the fit's GVD is 1.7 % above Malitson's formula's. The bands account for
  +2.1 % of that (the fundamental 1.3 %, the 2.2 µm group 0.8 %); the refit alone would
  be 0.3 % below. The 2.2 µm group alone changes β3 by 7 % at 2 µm, and near 2.2 µm
  the GVD swings through +700 fs²/mm. A finite-difference stencil a few nm wide
  resolves this. It is not an artefact of the fit. It is the dispersion of absorption
  measured in one sample:
  - Franta et al. do not assign the 2.2 µm bands;
  - their n swing, 1.2e-5 peak to peak, is about the precision of the absolute
    refractometry, whose wavelengths are 0.1-0.2 µm apart;
  - Malitson's three grades and Leviton's Corning 7980 agree with the fit to 1e-5 from
    1.5 to 2.44 µm. They neither confirm nor exclude bands of this strength in their
    glasses.

  A grade with less OH has weaker bands. Malitson's and Leviton's Sellmeier formulas
  have none (TC-841).
- **Evaluation noise:** the UDM's closed forms lose digits to cancellation below a
  term's characteristic energy (TC-321). n carries 1e-15 to 6e-15 of relative rounding
  noise from 0.8 to 3.5 µm, and dn/dω (a complex step, exact) 1e-13 to 6e-13. That is
  far below any physical effect, but it bounds a derivative differenced from them
  (TC-841).
- **Temperature:** Leviton & Frey's temperature-dependent formula (Corning 7980,
  30-300 K, 0.4-2.6 µm) is the model `leviton_2008/sellmeier`;
  `at_temperature` evaluates it at any temperature in that range. Their dn/dT (Table 6)
  is curated as datasets: 8.6e-6 per K at 800 nm and 295 K.

**In the infrared** the model's phonon bands are Franta's measurement of Lithosil Q2.
Philipp's older tables (Kramers-Kronig of reflectance) differ by up to 0.5 in n on the
reststrahlen bands (TC-324).

## The nonlinear index

The model (`kerr/SiO2`; TC-393) gives the long-pulse n2 and its Raman fraction f_R:

$$
n_{2,\mathrm{meas}}(\lambda, \tau) = n_2(\lambda)\left[(1 - f_R) + f_R\,\eta(\tau)\right],
\qquad n_2(\lambda) = N\left(1 + a E^2\right),
$$

- E is the photon energy (eV);
- eta(tau) is the share of the delayed Raman response a Gaussian pulse of duration tau
  sees, from Hollenbeck & Cantrell's 13-mode response;
- eta is 0.96 at 100 fs, 0.86 at 47 fs, and 1 for picosecond and longer pulses.

It is fitted to eleven sources' measured self coefficients, each with its own calibration
scale (a 10 % prior), from 350 nm to 3.5 µm, and to Stolen & Ippen's measured Raman gain
of fused quartz, which sets f_R (TC-696):

<!-- n2:start -->
| lambda (nm) | n2 (1e-20 m^2/W) | sigma | electronic |
|---|---|---|---|
| 400 | 3.05 | 0.18 | 2.39 |
| 532 | 2.78 | 0.13 | 2.17 |
| 800 | 2.59 | 0.11 | 2.02 |
| 1030 | 2.53 | 0.11 | 1.97 |
| 1550 | 2.47 | 0.11 | 1.93 |
| 2000 | 2.46 | 0.11 | 1.92 |
<!-- n2:end -->

The Raman fraction is f_R = 0.22 +- 0.03, set by Stolen & Ippen's gain, 1.5e-11 cm/W at
330 cm⁻¹ for a 526 nm pump, predicted as 2 omega_s n2 f_R Im H_R / c with Hollenbeck &
Cantrell's response (TC-695, TC-696). It agrees with Hellwarth et al.'s nuclear share of
n2 from Raman spectra, 0.21 +- 0.03, which is compared, not fitted. Stolen et al.'s
0.18, the prior until then, is left out: it was derived from these two measurements.
The electronic part is 0.78 of n2; a 50 fs pulse measures 0.97 of it. The Raman
response itself, with this f_R and two literature forms, is on {doc}`raman_response`.

- **Compared with the literature:**
  - Milam's much-used recommended values (2.74e-20 m^2/W at 1053 nm) are 8-11 % higher.
    His assessment rests largely on nanosecond measurements, which may include
    electrostriction.
  - Santran et al.'s ten grades are 37 % higher, beyond their stated absolute accuracy
    (10-28 %). Their pump-probe isolates the electronic part, which the larger f_R
    lowers.
  - The other nine sources lie within -13 % and +4 % of the fit. Milam & Weber's 1976
    interferometric value at 1064 nm (125 ps) is 4 % above it (TC-594).
- **Below 350 nm** two-photon absorption sets in (2 hbar omega above the 8.4 eV gap).
  DeSalvo et al.'s 266 nm value is twice this dispersion's extrapolation, and the model
  warns there.
- **Two-photon absorption, multiphoton and field ionisation, and laser damage:** their
  data are curated (TC-325 onwards). Their models are next.

## The nonlinear bundle

A propagation code in the $n_0 = 1$ convention of the nonlinear Schrödinger equation
takes the response $R(t) = (1 - f_R)\,\delta(t) + f_R\,h_R(t)$, with $\int h_R\,dt = 1$
and the total $n_2$: the triple $(n_2, f_R, h_R)$. `optmat.solid.bundle.silica_bundle`
gives it as one object whose parts cannot be mismatched (OEP 1, sections 4.7 and 6.9).
Every bundle holds the Kerr fit's **total** long-pulse $n_2(\lambda)$, the quantity the
data measure, and its shape's own $f_R$. Its electronic part is $(1 - f_R)\,n_2$, so a
scalar and a vector run agree in the long-pulse limit.

- `"default"`: Hollenbeck and Cantrell's shape with the Kerr fit's $f_R = 0.219$,
  `paired`, copolarised only.
- `"tensor"`: Lin and Agrawal's shape with its own $f_R = 0.245$ (TC-734), which has
  isotropic and anisotropic parts. Its electronic part is 3.3 % smaller than the
  default's, and its Raman part 12 % larger.
- One bundle for each silica response, under the response's name.

The silica response's `raman_coefficient` and `electronic_susceptibility` include
$n_0^2$; a code in the $n_0 = 1$ convention takes the bundle's parts instead.

Each bundle carries its gain check against Stolen and Ippen's measured gain,
$1.5 \times 10^{-11}$ cm/W ± 10 % at 330 cm⁻¹ for a 526 nm pump (`gain_check`). The Kerr
fit's $f_R$ was set by this gain through Hollenbeck's shape, so for another shape it is
a consistency check, not an independent validation:

| Bundle | $f_R$ | Gain at 330 cm⁻¹ (10⁻¹¹ cm/W) | Against the measurement |
|---|---:|---:|---:|
| `default` (Hollenbeck, paired) | 0.219 | 1.35 | −9.7 % |
| `tensor` (Lin and Agrawal, their $f_R$) | 0.245 | 1.38 | −7.9 % |
| `blow_1989` (Blow and Wood, their $f_R$) | 0.222 | 1.46 | −2.5 % |

```python
import numpy as np

from optmat.solid.bundle import silica_bundle

bundle = silica_bundle("tensor")
n2 = bundle.n2(1030e-9)  # m^2/W, the Kerr fit's total
fraction = bundle.raman_fraction
t = np.linspace(0.0, 1e-12, 2001)
h = bundle.h_r(t, "b")  # 1/s, the anisotropic part
assert abs(float(bundle.n2_electronic(1030e-9)) - (1 - fraction) * float(n2)) < 1e-30
assert abs(bundle.gain_check.deviation) < 0.1
```
