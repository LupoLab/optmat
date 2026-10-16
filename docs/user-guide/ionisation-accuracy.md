# Noble-gas ionisation: accuracy and limitations

Use these working uncertainty factors with the {doc}`/ionisation_rates` example. Their evidence is separated into {doc}`ionisation-static`, {doc}`ionisation-laser`, {doc}`ionisation-one-electron`, and {doc}`ionisation-measured`. The factors are judgements about model accuracy, not confidence intervals.

## Accuracy

<!-- ionisation-rate-accuracy:start -->
| Gas | Factor | Basis |
|---|---|---|
| He | 1.3 | Static: within 8 % of the all-electron rates from 0.07 to 0.15 au and 17 % to 0.3 au. 780 nm: 0.85 to 1.13 of Parker et al.'s fit from 1e14 to 1.4e15 W/cm2; 390 nm: 0.85 to 1.28 of it. One electron: 1.10 to 1.14 times the TDSE at 800 nm and 6 fs. At 390 nm with A2's 6 fs pulse Parker et al.'s rate integrated over the pulse is 0.97 to 1.00 of the two-electron calculation's ionisation plus excitation to n >= 5 (1.00 to 1.03 at its angular limit), which is 1.6 times one electron's at the same ionisation energy; haCC's helium is 0.96 (one ion channel) to 0.99 (five) of it, and 0.95 to 1.00 of MCTDHF at 800 nm (TC-660, TC-848, TC-1050, TC-1092). |
| Ne | 1.5 | Static: 0.82 to 0.96 of Jagau's all-electron rates from 0.14 to 0.24 au (TC-669). One electron at 400 nm: 1.14 times Zhao et al.'s calculation on geometric average (TC-831); no all-electron reference at finite frequency. |
| Ar | 1.5 | Static: 0.78 to 1.16 of Jagau's all-electron rates from 0.07 to 0.16 au (TC-669). One electron: 0.90 to 1.07 times the ionised fraction of the TDSE at 800 nm and 6 fs, 0.97 to 1.19 at 10 fs, 0.83 to 1.04 at 1030 nm, 0.59 to 1.36 at 400 nm (through resonances), 0.86 to 1.04 in Wallace et al.'s pulse (TC-848). The atom: 0.76 to 0.87 of one electron in haCC at 800 nm, about 0.87 in van der Hart's R-matrix rates at 390 nm (TC-842, TC-847). Wallace et al.'s yields are matched with the intensities scaled by 1.033, against 1.002 for the exact one-electron calculation (calibration 1.3 %; TC-836). Wahlstrand et al.'s absolute yields (42 fs) are 1.2 to 3.1 times lower than the model's, where the intensity's 13 % allows 1.9 (TC-660). |
| Kr | 2 | Static: one electron 1.08 to 1.22 of tRecX's and Bahl et al.'s one-electron rates; 1.04 of Santra et al.'s spin-orbit-resolved independent-electron sum (Herman-Skillman potential) at 0.10 au (TC-821, TC-845, TC-846). One electron: 1.16 (6 fs) and 1.20 (Wallace et al.'s pulse) times the TDSE at 800 nm, the static normalisation included. Wallace et al.'s yields are matched with the intensities scaled by 1.023 (calibration 1.5 %); Wahlstrand et al.'s are 1.4 to 3.8 times lower than the model's (TC-660, TC-848). |
| Xe | 3 | Static: one electron 0.84 to 1.05 of tRecX's and Bahl et al.'s one-electron rates (TC-821, TC-845). One electron: 0.94 times the TDSE at 800 nm and 6 fs, 0.96 to 1.00 in Wallace et al.'s pulse (TC-849). The atom's 1.5 electrons rest on Wallace et al.'s yields alone (scale 1.026, calibration 2.5 %). Against it, within the absolute measurements at 800 nm, where one calibration serves all gases, xenon's intensity scale is 0.86 (Nambu et al.) and 0.81 (Nie et al.) of argon's, the model's xenon high against its argon by a factor of 2 to 3.5 in the rate (TC-819); at 400 nm Nambu et al.'s is 1.09. Wahlstrand et al.'s yields are 2.9 to 6.8 times lower than the model's. The structure coefficients of the sources differ by a factor of 3 in C^2 for xenon (TC-657, TC-660, TC-850). Hartmann et al.'s all-electron R-matrix (RMT) calculation in a 5.3 fs pulse at 759 nm depletes the ground state like 0.53 to 0.64 of one m = 0 electron at its ionisation energy, and the default 1.5 to 2.1 times as much as RMT; against the same paper's measured neutral abundance the default's depletion is 0.96 to 1.13 up to 1.6e14 W/cm2 (TC-1093, TC-1094). |
<!-- ionisation-rate-accuracy:end -->

The factors are judgements from the comparisons of the preceding pages, not statistical intervals.
They apply between $10^{-6}$ of saturation and saturation. The yield of a pulse goes
as the fifth to eighth power of the intensity near 1e14 W/cm², so a factor of 1.5 in
the rate is 5 to 8 % in the intensity.

## Open questions

- **Short wavelengths and two photons.** Above half the ionisation energy, a
  vacuum-ultraviolet field ionises by absorbing two of its own photons, which the
  strong-field models (built at a carrier, 350 nm and longer) cannot represent; fed a
  total field, they see a vacuum-ultraviolet component only as a ripple on $|E|$. The
  two-photon model is {mod}`optmat.gas.two_photon` ({doc}`/two_photon`; OEP 6). Applying
  both it and a strong-field rate to the same spectral component counts that component
  twice; the band split is the propagation code's.
- **The two measured yields** differ by about 3 in the rate ({doc}`ionisation-measured`). The default
  follows Wallace et al.'s in the meantime (TC-814).
- **Krypton's effective number of electrons** (1.35) rests on Wallace et al.'s yields
  through the exact one-electron calculation, which want about 1.4 (TC-846, TC-851).
  The one all-electron calculation, Pathak et al.'s (2022) time-dependent CASSCF of the
  eight 4s4p electrons at 800 nm and 2e14 W/cm2, gives 0.91 of one electron, without
  spin-orbit coupling and in an assumed pulse envelope (TC-860). Not reconciled; left
  as it is (the maintainer, 2026-10-06).
- **Xenon against argon.** Three absolute measurements put the static-factor default's
  xenon a factor of 2 to 3 high against its argon (TC-819); see also *Xenon's effective
  number of electrons* below.
- **The Coulomb factor at large Keldysh parameter.** At 400 nm the static factor was 3
  to 5 times above Nambu et al.'s measured ionisation of Ar, Kr and Xe (TC-817); the
  frequency-dependent factor, the default since TC-848, removes most of it. Above
  $\gamma = 2$, in the multiphoton regime, it is low for hydrogen (TC-833).
- **The multiphoton regime.** At 248.6 nm, below the models' range, the default's
  rate is 3.4 and 4.3 times what helium's and argon's measured saturation intensities
  imply, and krypton's 0.55 (TC-853); for hydrogen above gamma = 2 it is low. Part of
  it is the integral over the channels: with the sum, argon's all-electron
  three-photon cross section at 193 nm is matched within 26 % (TC-856). The default
  keeps the integral; the sum is the option for the multiphoton regime.
- **Lai et al.'s yields** (Ne, Ar, Kr, Xe from 0.4 to 3.9 um, linear and circular;
  TC-857) are curated but not yet compared: their intensities were in some sets
  rescaled by up to 20 % to fit PPT, so they test the shape of the yield (its slope,
  its saturation, circular over linear), not the rate's magnitude.
- **Xenon's effective number of electrons** (1.5) rests on Wallace et al.'s yields
  alone, through the one-electron calculation in their pulse (TC-849, TC-850). Within
  the other absolute measurements at 800 nm, where one calibration serves all gases,
  the model's xenon is high against its argon (intensity scales 0.86 and 0.81 of
  argon's, a factor of 2 to 3.5 in the rate); at 400 nm Nambu et al.'s is 1.09 (TC-819).
  No all-electron calculation of xenon decides it.
- **The Coulomb correction in circular polarisation** is borrowed from linear
  polarisation (TC-820).
- **The earlier static factor.** Before TC-848, for one electron the model was 1.3 to 1.7 times above
  the calculations at Keldysh parameters 1.2 to 2.3, and none of the then-compared analytic Coulomb factors
  followed them (TC-812). The empirical frequency-dependent factor addresses this
  excess ({doc}`ionisation-one-electron`); its extrapolation remains uncertain. Wallace et al.'s argon tests the rate where the yield bends
  towards saturation, near 3e14 to 5e14 W/cm², where the model is only 1.07 to 1.12
  times the calculation; at lower intensities, where it is up to 1.7 times, the
  measurements of the absolute yield are Wahlstrand et al.'s alone.
- **The effective number of electrons at low field.** The all-electron static rates
  of argon rise above the model below 0.07 au (1.7e14 W/cm²), and the measured yields
  of Wahlstrand et al. lie below it there. No constant satisfies both.
- **Krypton and xenon** have one-electron static references (Bahl et al., within 0.88
  to 1.22) and, for krypton, one spin-orbit-resolved independent-electron sum
  (Santra et al.: 0.77 with argon's one electron, 1.04 with krypton's 1.35), but no
  all-electron static rates; the later krypton laser calculation and its unresolved
  interpretation are noted above. Their structure coefficients differ between
  sources by a factor of up to 3 in $C_l^2$ (xenon) (TC-821).
- **400 nm and shorter.** The earlier static-factor model at 390 nm was 27 to 42 % above the two-electron
  rates of helium, and 1.6 to 2.9 times above the all-electron R-matrix rates of argon
  (7 to 12 times where the rate is multiphoton), though its one-electron form matches a
  one-electron calculation of argon (TC-701). The frequency-dependent factor removes
  most of this excess: the default is 0.82 to 1.30 times the argon R-matrix rates, but
  still 2.6 to 3.4 times where the rate is multiphoton ({doc}`ionisation-laser`).
  Whether the argon atom's rate is below one active electron's, and whether at 800 nm
  too, is open. Resonances and channel closings, which the rate does not contain,
  matter most there.
- **Molecules** are covered on {doc}`/molecular_ionisation_rates`.
- **1030 nm and longer.** One electron of argon at 1030 nm (10 fs) is within 1.35 of
  the older static-factor model (TC-826), and the same form is 1.09 times Zhao et al.'s hydrogen at 1200
  nm (TC-827). No measurement there is curated.
- **Hydrogen.** Popruzhenko et al. validated their Coulomb factor on hydrogen at
  400 nm. Against Zhao et al.'s one-electron hydrogen at 600 to 1200 nm the static
  factor, with Tong and Lin's static correction (alpha = 6), is 1.09 to 1.53 times
  high (TC-827).
- **The timing** of the ionisation within the cycle at $\gamma \gtrsim 1$.

## Comparison with a recorded Luna.jl snapshot

This comparison uses Luna.jl at commit `de9db7c`, recorded in
`tests/data/legacy_ionisation_rates.json` (TC-815). It describes that snapshot.
The rates were evaluated on the same field grid at 800 nm; the first table's yields
use a 30 fs Gaussian intensity envelope. The snapshot uses PPT in the form of Ilkov
et al. and Bergé et al., with $(1 + \gamma^2)^{-n^*}$ in the Coulomb factor, a Hartree
coefficient with $l^* = n^* - 1$, a sum over $m$, and two electrons per orbital.
For He, Ne and Ar it also includes a Stark shift and $\exp(-2\alpha_+F)$.
It has no Tong–Lin finite-field correction and refuses fields above twice the
barrier-suppression field; blank entries mark that limit.

<!-- ionisation-rate-luna:start -->
| Gas | I (W/cm2) | gamma | optmat / Luna | 30 fs yield: Luna | optmat |
| --- | --- | --- | --- | --- | --- |
| He | 5.00e+13 | 2.03 | 1.62 | 2.8e-12 | 4.7e-12 |
| He | 1.00e+14 | 1.43 | 1.24 | 8.5e-09 | 1.1e-08 |
| He | 2.00e+14 | 1.01 | 0.99 | 7.4e-06 | 7.6e-06 |
| He | 5.00e+14 | 0.64 | 0.79 | 0.0073 | 0.006 |
| Ne | 5.00e+13 | 1.90 | 0.50 | 3e-09 | 1.6e-09 |
| Ne | 1.00e+14 | 1.34 | 0.38 | 2.8e-06 | 1.1e-06 |
| Ne | 2.00e+14 | 0.95 | 0.30 | 0.00084 | 0.00026 |
| Ne | 5.00e+14 | 0.60 | 0.24 | 0.22 | 0.061 |
| Ar | 5.00e+13 | 1.62 | 1.13 | 9.5e-05 | 0.00011 |
| Ar | 1.00e+14 | 1.15 | 0.97 | 0.0081 | 0.0081 |
| Ar | 2.00e+14 | 0.81 | 0.92 | 0.24 | 0.22 |
| Ar | 5.00e+14 | 0.51 | 1.03 | 1 | 1 |
| Kr | 5.00e+13 | 1.53 | 0.67 | 0.0039 | 0.0029 |
| Kr | 1.00e+14 | 1.08 | 0.41 | 0.21 | 0.1 |
| Kr | 2.00e+14 | 0.77 | 0.25 | 1 | 0.84 |
| Kr | 5.00e+14 | 0.48 | 0.11 | 1 | 1 |
| Xe | 5.00e+13 | 1.42 | 0.64 | 0.09 | 0.068 |
| Xe | 1.00e+14 | 1.01 | 0.36 | 0.94 | 0.71 |
| Xe | 2.00e+14 | 0.71 | 0.19 | 1 | 1 |
| Xe | 5.00e+14 | 0.45 |  |  | 1 |
<!-- ionisation-rate-luna:end -->

The rate ratio depends on gas and intensity. In this grid it is 0.92–1.13 for Ar,
0.79–1.62 for He, and below one for Ne, Kr and Xe. These are differences between
models, not uncertainty estimates.

The following tests use each rate in the same yield calculations as the preceding
validation chapters (TC-815, TC-829). Recorded rates are interpolated on the
fixture's field grid. Each row states its observable: an intensity scale or a rate
ratio. Parker et al.'s helium reference is at 780 nm; the main comparison is at
800 nm, with optmat's 780 nm result also shown.

<!-- ionisation-rate-luna-tests:start -->
| Test | Gas | Luna | optmat |
| --- | --- | --- | --- |
| Wallace et al. 2016: intensity factor (calibration 1.3 %) | Ar | 1.05 | 1.03 |
| Wallace et al. 2016: intensity factor (calibration 1.5 %) | Kr | 0.95 | 1.02 |
| Wallace et al. 2016: intensity factor (calibration 2.5 %) | Xe | 0.90 | 1.03 |
| `nambu_2026/fig6_Ar_LP`: intensity factor | Ar | 1.08 | 1.06 |
| `nambu_2026/fig6_Kr_LP`: intensity factor | Kr | 0.85 | 0.97 |
| `nambu_2026/fig6_Xe_LP`: intensity factor | Xe | 0.88 | 0.92 |
| `nie_2022/fig3_He`: intensity factor | He | 0.90 | 0.91 |
| `nie_2022/fig4_Ar_LP`: intensity factor | Ar | 0.98 | 0.98 |
| `nie_2022/fig4_Kr_LP`: intensity factor | Kr | 0.80 | 0.89 |
| `nie_2022/fig4_Xe_LP`: intensity factor | Xe | 0.75 | 0.80 |
| Wahlstrand et al. 2018: rate over measured | Ar | 0.97 to 3.10 | 1.20 to 3.06 |
| Wahlstrand et al. 2018: rate over measured | Kr | 1.44 to 7.31 | 1.36 to 3.80 |
| Wahlstrand et al. 2018: rate over measured | Xe | 1.81 to 9.37 | 2.90 to 6.75 |
| Parker et al. 2007 (780 nm, two electrons, to 1e15 W/cm2): rate over theirs, the codes at 800 nm | He | 0.59 to 1.44 | 0.76 to 1.05 (at 780 nm: 0.85 to 1.06) |
<!-- ionisation-rate-luna-tests:end -->

The intensity factors depend on the measurement as well as the rate. Across the ten
intensity-scale fits the root-mean-square factor is 1.16 for this Luna snapshot and
1.10 for optmat's default. The disagreement between Wallace, Wahlstrand, Nambu and
Nie remains part of the uncertainty assessment above; agreement with one series
does not settle the others. The effective electron counts of Kr and Xe were selected
using Wallace et al.'s yields, so that comparison is not an independent validation.
