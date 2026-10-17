# Noble-gas ionisation: all-electron laser calculations

These comparisons test the atom in a laser field. Read {doc}`ionisation-static` for the static normalisation and {doc}`ionisation-theory` for the distinction between cycle-averaged and instantaneous rates. The separate {doc}`ionisation-one-electron` chapter isolates the frequency dependence using model potentials.

## Helium in a laser field

Parker et al.'s (2007) fits of the rates of their two-electron calculations,
at 780 and 390 nm. The rate is that of the atom. The columns are the default's rate over
theirs, and the same with the other Coulomb factors (the static factor's column was
labelled `default` until TC-1050).

<!-- ionisation-rate-helium:start -->
| Wavelength | I (W/cm2) | gamma | Parker (au) | `default` | with the static factor | with Ilkov's factor | with Popruzhenko's |
|---|---|---|---|---|---|---|---|
| 780 nm | 1.00e+14 | 1.47 | 3.53e-11 | 0.87 | 1.00 | 0.42 | 0.34 |
| 780 nm | 2.00e+14 | 1.04 | 2.11e-08 | 0.86 | 0.95 | 0.55 | 0.41 |
| 780 nm | 3.00e+14 | 0.85 | 4.67e-07 | 0.89 | 0.98 | 0.65 | 0.47 |
| 780 nm | 4.00e+14 | 0.74 | 3.17e-06 | 0.93 | 1.01 | 0.73 | 0.53 |
| 780 nm | 6.00e+14 | 0.60 | 3.27e-05 | 0.99 | 1.07 | 0.85 | 0.62 |
| 780 nm | 8.00e+14 | 0.52 | 1.34e-04 | 1.04 | 1.11 | 0.93 | 0.68 |
| 780 nm | 1.00e+15 | 0.47 | 3.55e-04 | 1.08 | 1.14 | 0.98 | 0.73 |
| 780 nm | 1.40e+15 | 0.39 | 1.26e-03 | 1.13 | 1.19 | 1.07 | 0.81 |
| 390 nm | 1.00e+14 | 2.94 | 4.95e-09 | 0.85 | 1.31 | 0.24 | 0.24 |
| 390 nm | 2.00e+14 | 2.08 | 2.72e-07 | 1.01 | 1.43 | 0.41 | 0.36 |
| 390 nm | 3.00e+14 | 1.70 | 2.54e-06 | 1.00 | 1.36 | 0.49 | 0.41 |
| 390 nm | 4.00e+14 | 1.47 | 1.11e-05 | 1.00 | 1.31 | 0.56 | 0.44 |
| 390 nm | 6.00e+14 | 1.20 | 7.22e-05 | 1.01 | 1.28 | 0.66 | 0.50 |
| 390 nm | 8.00e+14 | 1.04 | 2.35e-04 | 1.03 | 1.27 | 0.74 | 0.55 |
| 390 nm | 1.00e+15 | 0.93 | 5.37e-04 | 1.06 | 1.28 | 0.81 | 0.59 |
| 390 nm | 1.40e+15 | 0.79 | 1.62e-03 | 1.11 | 1.31 | 0.92 | 0.67 |
<!-- ionisation-rate-helium:end -->

**Two electrons at 390 nm** (TC-1050, TC-1092): helium with both electrons active in tRecX,
A2's 6 fs pulse at 390 nm, against haCC's helium (one ion channel) and one electron at
the same ionisation energy, all read as the norm-loss 1000 au after the pulse
({func}`~optmat.gas.ionisation_validation.helium_two_electron_ratios`), and Parker et
al.'s rate integrated over the same pulse
({func}`~optmat.gas.ionisation_validation.parker_pulse_yield`):

<!-- ionisation-rate-helium-two-electron:start -->
| Calculation | I (W/cm2) | Norm-loss | One electron, same Ip | Over one electron | Parker et al. over the pulse | Over Parker et al. |
|---|---|---|---|---|---|---|
| `twoelectron_normloss_gauss6_390nm_He` | 5.69e+14 | 0.00733 | 0.004289 | 1.71 | 0.007 | 1.05 |
| `twoelectron_normloss_gauss6_390nm_He` | 3.68e+14 | 0.0008623 | 0.0004807 | 1.79 | 0.0008491 | 1.02 |
| `twoelectron_normloss_gauss6_390nm_He_m5` | 5.69e+14 | 0.007207 | 0.004186 | 1.72 | 0.007 | 1.03 |
| `hacc_normloss_gauss6_390nm_He` | 5.69e+14 | 0.005595 | 0.004289 | 1.30 | 0.007 | 0.80 |
| `hacc_normloss_gauss6_390nm_He` | 3.68e+14 | 0.0006777 | 0.0004807 | 1.41 | 0.0008491 | 0.80 |
| `twoelectron_normloss_gauss6_390nm_He_box60` | 5.69e+14 | 0.006901 |  |  | 0.007 | 0.99 |
<!-- ionisation-rate-helium-two-electron:end -->

- **These norm-losses are not one measure.** The two-electron norm-loss falls as the
  box grows (5.9 % from BOX 40 to 60), because each larger box keeps another Rydberg
  shell bound; haCC's bound states are converged at 40 bohr, so its norm-loss is the
  ionisation alone. The ratios above therefore set ionisation plus Rydberg excitation
  against ionisation, and the two-electron basis lies about 3 % above its angular limit
  (TC-1092). They are kept as the record of these runs; the comparison is below.

**Matched measures** (TC-1092): the same atom with both electrons active, haCC with one
ion channel (haCC(1)) and with the five n ≤ 2 channels (haCC(5)), and one electron, in
the same pulse, read as the population above the bound shells n ≤ 4 (ionisation plus
excitation to n ≥ 5) or n ≤ 5, 1000 au after the pulse: measures converged in the box
(`travers_2026_pyroraptor`). At 800 nm, haCC against Sato et al.'s MCTDHF in their
pulse and with their measure
({func}`~optmat.gas.ionisation_validation.helium_matched_ratios`):

<!-- ionisation-rate-helium-matched:start -->
| Calculation | Against | I (W/cm2) | Yield | Against's | Ratio | ln(1 - Y) ratio |
|---|---|---|---|---|---|---|
| `twoelectron_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/sae_above4_gauss6_390nm_He` | 3.68e+14 | 0.0008464 | 0.0005119 | 1.653 | 1.654 |
| `twoelectron_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/sae_above4_gauss6_390nm_He` | 5.69e+14 | 0.007244 | 0.004504 | 1.608 | 1.611 |
| `hacc1_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 3.68e+14 | 0.0008092 | 0.0008464 | 0.956 | 0.956 |
| `hacc1_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 5.69e+14 | 0.006985 | 0.007244 | 0.964 | 0.964 |
| `hacc5_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 3.68e+14 | 0.0008277 | 0.0008464 | 0.978 | 0.978 |
| `hacc5_above4_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 5.69e+14 | 0.00715 | 0.007244 | 0.987 | 0.987 |
| `hacc1_above5_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He` | 3.68e+14 | 0.0007113 | 0.0007396 | 0.962 | 0.962 |
| `hacc1_above5_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He` | 5.69e+14 | 0.006513 | 0.006742 | 0.966 | 0.966 |
| `hacc5_above5_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He` | 3.68e+14 | 0.0007242 | 0.0007396 | 0.979 | 0.979 |
| `hacc5_above5_gauss6_390nm_He` | `travers_2026_pyroraptor/twoelectron_above5_gauss6_390nm_He` | 5.69e+14 | 0.006672 | 0.006742 | 0.990 | 0.990 |
| `hacc1_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 4.00e+14 | 0.0002961 | 0.0003029 | 0.978 | 0.978 |
| `hacc1_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 8.00e+14 | 0.01545 | 0.01589 | 0.972 | 0.972 |
| `hacc1_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 1.40e+15 | 0.1536 | 0.1608 | 0.955 | 0.951 |
| `hacc5_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 4.00e+14 | 0.0003023 | 0.0003029 | 0.998 | 0.998 |
| `hacc5_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 8.00e+14 | 0.01583 | 0.01589 | 0.996 | 0.996 |
| `hacc5_outside20_sin2_800nm_He` | `sato_2016/fig6_casscf_He_800nm` | 1.40e+15 | 0.1578 | 0.1608 | 0.981 | 0.979 |
| Parker et al.'s rate over the pulse | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 3.68e+14 | 0.0008491 | 0.0008464 | 1.003 | 1.003 |
| Parker et al.'s rate over the pulse | `travers_2026_pyroraptor/twoelectron_above4_gauss6_390nm_He` | 5.69e+14 | 0.007 | 0.007244 | 0.966 | 0.966 |
<!-- ionisation-rate-helium-matched:end -->

- **The atom ionises like 1.6 electrons** of one electron at its own ionisation energy:
  1.65 and 1.61, or 1.60 and 1.56 at the two-electron calculation's angular limit (3.1
  % lower). Sato et al.'s converged helium at 800 nm gives 1.71 to 1.53, falling with
  the intensity (TC-860). A single number for helium is uncertain by about 10 %.
- **Parker et al.'s rate is confirmed:** integrated over the pulse it is 1.00 and 0.97
  of the two-electron yield (1.03 and 1.00 at the angular limit). The default's
  $N_\text{eff}$ of 1.47, set from their static rates, stands.
- **haCC's helium is within a few per cent of the exact atom:** haCC(5) is 0.98 and
  0.99 of the two-electron calculation at 390 nm (above the angular limit by about the
  calculation's 3 %) and 0.98 to 1.00 of MCTDHF at 800 nm; haCC(1) 0.96 at 390 nm and
  0.95 to 0.98 at 800 nm. The 0.76 to 0.79 above, and the factor 1.25 between haCC's
  norm-loss and Parker et al.'s rate (TC-848), came from the measure, not from haCC.
  Whether haCC's argon with the Ar_qz data is as close is not known.
- **Convergence:** the populations above n ≤ 4 change by 0.03 % from R0 = 78 to 97 bohr
  (two electrons) and by 0.3 % among haCC's numerical variations; the two-electron
  calculation's angular basis (3.1 %) is the largest uncertainty.

## All-electron rates at 390 nm

Van der Hart's R-matrix calculations describe the whole atom, with one state of the
ion: R-matrix Floquet rates of He, Ne and Ar (2006; read from the image of its Fig. 6)
and time-dependent R-matrix rates of Ar (2007; its Table I, and the rest of its Fig. 4
read from the vector drawing). The figures are read by data rule 8: the Ar curve read
from the 2006 image matches the 2007 vector copy of it within 4 %, and the 2007 figure
matches its own Table I within 6 % (TC-700, TC-701). The references resolve the atoms'
resonances, which no model has; "Parker et al." is the fit of helium's two-electron
calculation over the same reference.

<!-- ionisation-rate-rmatrix:start -->
| Gas | Reference | I (W/cm2) | gamma | Reference (1/s) | `default` | `one_electron_static` | `one_electron` (the default's factor) | Parker et al. |
|---|---|---|---|---|---|---|---|---|
| He | `vanderhart_2006/fig6_floquet_He` | 1.10e+14 | 2.81 | 3.33e+08 | 0.94 | 0.97 | 0.64 | 1.05 |
| He | `vanderhart_2006/fig6_floquet_He` | 1.50e+14 | 2.40 | 6.09e+09 | 0.35 | 0.35 | 0.24 | 0.35 |
| He | `vanderhart_2006/fig6_floquet_He` | 2.00e+14 | 2.08 | 1.18e+10 | 0.97 | 0.93 | 0.66 | 0.96 |
| He | `vanderhart_2006/fig6_floquet_He` | 2.40e+14 | 1.90 | 2.51e+10 | 1.26 | 1.19 | 0.86 | 1.25 |
| Ne | `vanderhart_2006/fig6_floquet_Ne` | 1.10e+14 | 2.63 | 4.72e+09 | 1.84 | 3.22 | 1.90 |  |
| Ne | `vanderhart_2006/fig6_floquet_Ne` | 1.50e+14 | 2.25 | 5.20e+10 | 0.84 | 1.40 | 0.87 |  |
| Ne | `vanderhart_2006/fig6_floquet_Ne` | 2.00e+14 | 1.95 | 1.23e+11 | 1.45 | 2.31 | 1.50 |  |
| Ar | `vanderhart_2006/fig6_floquet_Ar` | 6.50e+13 | 2.92 | 2.46e+11 | 3.37 | 10.24 | 3.59 |  |
| Ar | `vanderhart_2006/fig6_floquet_Ar` | 1.00e+14 | 2.35 | 9.28e+12 | 0.43 | 1.15 | 0.46 |  |
| Ar | `vanderhart_2006/fig6_floquet_Ar` | 1.50e+14 | 1.92 | 1.22e+13 | 1.30 | 3.05 | 1.36 |  |
| Ar | `vanderhart_2006/fig6_floquet_Ar` | 1.90e+14 | 1.71 | 2.73e+13 | 1.20 | 2.64 | 1.25 |  |
| Ar | `vanderhart_2007/table1_one_state_Ar` | 4.80e+13 | 3.40 | 9.79e+10 | 2.60 | 8.76 | 2.78 |  |
| Ar | `vanderhart_2007/table1_one_state_Ar` | 7.50e+13 | 2.72 | 4.32e+11 | 3.27 | 9.50 | 3.48 |  |
| Ar | `vanderhart_2007/table1_one_state_Ar` | 1.50e+14 | 1.92 | 1.21e+13 | 1.30 | 3.06 | 1.36 |  |
| Ar | `vanderhart_2007/fig4_td_rmatrix_Ar` | 2.50e+13 | 4.71 | 5.53e+09 | 3.18 | 13.51 | 3.44 |  |
| Ar | `vanderhart_2007/fig4_td_rmatrix_Ar` | 1.25e+14 | 2.11 | 9.28e+12 | 0.93 | 2.31 | 0.98 |  |
| Ar | `vanderhart_2007/fig4_td_rmatrix_Ar` | 1.75e+14 | 1.78 | 2.02e+13 | 1.26 | 2.84 | 1.32 |  |
| Ar | `vanderhart_2007/fig4_td_rmatrix_Ar` | 2.50e+14 | 1.49 | 8.86e+13 | 0.82 | 1.67 | 0.85 |  |
| Ar | `vanderhart_2007/fig4_td_rmatrix_Ar` | 3.00e+14 | 1.36 | 1.03e+14 | 1.13 | 2.22 | 1.18 |  |
<!-- ionisation-rate-rmatrix:end -->

- **Helium:** away from the resonance at 1.5e14 W/cm², Parker's fit is 0.96–1.25
  times the Floquet rate and the current default is 0.94–1.26 times it. The earlier
  static-factor default was 1.4–1.8 times above it (TC-701).
- **Argon:** at 1.25e14–3e14 W/cm² the current default is 0.82–1.30 times the
  listed time-dependent and Floquet rates. At the lower-intensity multiphoton points
  it remains 2.6–3.4 times high. The older static-factor comparison was 1.6–2.9 times
  high at the higher intensities and 7–12 times at the lower ones (TC-701). The
  frequency-dependent correction removes much of that excess. Its one-electron
  rate is 1.17 times the Floquet rates and 1.07 times the time-dependent rates on
  geometric average below $\gamma = 2.2$ (0.85–1.37). The inferred atom-to-electron
  ratio, about 0.87, is consistent with the interim haCC comparison at 800 nm
  (0.76–0.87; TC-830, TC-842).
- **Neon:** the current default is 0.84–1.84 times the three listed rates. The older
  static-factor comparison gave an excess of about 1.4–3.

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_ionisation_rmatrix_001.*
:target: ../auto_gallery/ionisation/plot_ionisation_rmatrix.html
:alt: The noble gases' default rate against van der Hart's all-electron rates at 390 nm

Van der Hart's R-matrix rates of He, Ne and Ar at 390 nm with the default's; the
default over them against the Keldysh parameter, with Parker et al.'s helium over the
Floquet rate.
```
