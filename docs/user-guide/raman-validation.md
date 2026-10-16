# Raman validation and limitations

Read {doc}`/raman` and {doc}`/raman_response` first. The comparisons below test line
strengths, widths and assembled kernels against measured gains and refractive-index
responses. They retain the distinction between source uncertainty and model error.
The source catalogue is {doc}`/the-data/sources`.

## How well the lines give the measured gains

The Raman gain tests the strengths and the widths together, and every value here was
measured independently of both. The table gives the model's gain at 295 K for the
line and density of each measurement:

<!-- raman-gain:start -->
| Gas | Line | Pump (nm) | Density (amagat) | Measurement | Measured (cm/GW) | Model (cm/GW) | Model / measured - 1 |
|---|---|---|---|---|---|---|---|
| H2 | Q(1) | 532 | 60 | ottusch_1988/gain_H2 | 2.64 +- 0.20 | 3.36 | +27% |
| H2 | Q(1) | 532 | 20 | bischel_1986b/gain_H2 | 2.50 +- 0.40 | 3.31 | +33% |
| H2 | Q(1) | 477 | 20 | bischel_1986b/gain_H2 | 3.50 +- 0.30 | 3.90 | +11% |
| H2 | Q(1) | 308 | 20 | bischel_1986b/gain_H2 | 6.60 +- 0.80 | 7.95 | +21% |
| D2 | Q(2) | 532 | 60 | ottusch_1988/gain_D2 | 0.45 +- 0.05 | 0.62 | +39% |
| CH4 | Q | 532 | 115 | ottusch_1988/gain_CH4 | 1.26 +- 0.10 | 1.26 | +0% |
<!-- raman-gain:end -->

- **H2 and D2: the strengths.** These are now ab initio (Hunt et al.) and dispersed to
  the pump. They agree with Rychlewski's frequency-dependent elements to better than
  1 %.
- **H2 and D2: the gains.** The model's gains are 11-33 % above the H2 measurements and
  39 % above D2's. Before, a single 500 nm element was held at every pump, and the
  model ranged from 10 % low at 308 nm to 23 % high at 532 nm. The error that varied
  with the pump is gone, but a constant excess remains.
  - The excess is therefore not in the strengths. It lies in the widths or in the
    measured gains.
  - Ottusch & Rockwell's and Bischel & Dyer's own calculated gains are 8 % and 18 %
    above their measurements at 532 nm.
  - Nothing is tuned to close this gap (TC-542).
- **CH4:** with the strength from Schrötter & Klöckner's Σ, dispersed to the pump, the
  model agrees at 532 nm. With Murphy et al.'s strength it was 7 % low.

## The rotational nonlinear index

A pulse long against the rotational periods sees the time integral of the rotational
response, and its index changes by n2,rot I with n2,rot = N Σ A_J/Ω_J/(2ε₀c)
({func}`~optmat.gas.raman.rotational_nonlinear_index`). The response is slow, so a
probe sees the same shift as the pulse itself; for the instantaneous electronic
response the probe's shift is twice the self shift. Against the published values at
1 atm:

<!-- raman-rotational:start -->
| Gas | Published | T (K) | n2,rot published (1e-20 cm2/W) | Model | Model / published - 1 |
|---|---|---|---|---|---|
| N2 | zahedpour_2015/n2_rot_N2 | 293.15 | 24 | 23.9 | -0% |
| N2 | zahedpour_2019/n2_rot_N2 | 293.15 | 23 | 23.9 | +4% |
| O2 | zahedpour_2015/n2_rot_O2 | 293.15 | 54 | 55.2 | +2% |
| O2 | zahedpour_2019/n2_rot_O2 | 293.15 | 53.3 | 55.2 | +4% |
| H2 | wahlstrand_2015/n2_rot_H2 | 290.9 | 2.7 | 3.1 | +15% |
| D2 | wahlstrand_2015/n2_rot_D2 | 290.9 | 3.9 | 4.0 | +2% |
<!-- raman-rotational:end -->

- **The factor of two of TC-312 is settled (TC-493).** The lines reproduce the
  published n2,rot of N2 and O2 to 0-4 % and of D2 to 2 %, taken as self coefficients.
  These lines use the same anisotropies as the measurements' calibrations. Zahedpour
  et al.'s (2019) printed response function integrates to twice their n2,rot (TC-312),
  so the factor of two lies in that normalisation, not in the n2,rot.
- **H2:** 15 % high. optmat uses the ab initio anisotropy, which is 4.5 % above
  Wahlstrand et al.'s measured one (squared, 9 %), and depends on J.
- **The classical limit:** the coherent (ΔJ = ±2) response is 3/4 of the classical,
  isothermal orientational response Δα²/(22.5 ε₀k_BT). The remaining quarter is the
  redistribution of the M sublevels, which collisions establish over nanoseconds and a
  short pulse does not.


## Line-model open questions

- **Level energies:** modern level energies for H2 and D2 (Komasa et al.), and the
  rotational constants of N2O. Both are requested.
- **The H2 and D2 gains:** the 11-39 % excess over the measured gains (see above).
- **Stopgaps:**
  - the Q-branch widths of CO2 and N2O;
  - N2O's depolarisation ratios;
  - one Σ for SF6, whose curated values span 20 %.
- **Line mixing:** the collapse of the Q branch above a few atmospheres. Rosasco et
  al.'s mixing coefficients are recorded, but N2's fail the detailed-balance sum rule
  (TC-527).
- **The dispersion of the anisotropies:** Hunt et al.'s static anisotropies are 3-4 %
  below the 632.8 nm values used.

## Interpreting older parameterisations

A width assigned to an unresolved Q branch is not an individual line's collisional
width. For N2, a 6 ps vibrational dephasing time describes the room-temperature
branch as a whole; the resolved model gives each J its own width (TC-491).
Wahlstrand et al.'s H2 and D2 rotational constants, 58.9 and 29.3 cm⁻¹, are effective
fit values. Their vibrational derivatives (1.3 and 1.4e-16 cm², ±15–20 %) give
amplitudes 13 % and 31 % larger than the curated static Q(1) and Q(2) elements;
the Raman-scattering value also quoted for D2, 1.22, agrees with Murphy et al.
and Hunt et al. A single representative line cannot describe Q-branch beating.
Murphy et al.'s CH4 and SF6 strengths give amplitudes 10 % and 11 % below those
from Schrötter & Klöckner's cross sections (TC-518, TC-542).

## Response validation

**Gains.** The gains from the imaginary part of the response at each line (its
neighbours and its anisotropic share included) against measurements that entered none
of its parameters. They are those of the line gains above, 2-3 % higher for H2 and D2, by the
anisotropic share of Q(1) and Q(2):

<!-- response-gain:start -->
| Gas | Line | Measurement | Measured (cm/GW) | Response (cm/GW) | Response / measured - 1 |
|---|---|---|---|---|---|
| H2 | Q(1), 532 nm, 60 amagat | ottusch_1988/gain_H2 | 2.64 +- 0.2 | 3.44 | +30% |
| H2 | Q(1), 532 nm, 20 amagat | bischel_1986b/gain_H2 | 2.5 +- 0.4 | 3.39 | +36% |
| H2 | Q(1), 477 nm, 20 amagat | bischel_1986b/gain_H2 | 3.5 +- 0.3 | 4 | +14% |
| H2 | Q(1), 308 nm, 20 amagat | bischel_1986b/gain_H2 | 6.6 +- 0.8 | 8.15 | +23% |
| D2 | Q(2), 532 nm, 60 amagat | ottusch_1988/gain_D2 | 0.45 +- 0.05 | 0.635 | +41% |
| CH4 | Q, 532 nm, 115 amagat | ottusch_1988/gain_CH4 | 1.26 +- 0.1 | 1.26 | +0% |
<!-- response-gain:end -->

The H2 and D2 excess is TC-542's: it is in neither the strengths nor the response's
assembly, and nothing is tuned to close it.

**The rotational revivals.** Chen et al. fitted an exponential to the declining
revival peaks of O2 at 5.1 atm; the same fit to the response's revivals, with
Hölzer et al.'s per-N widths:

<!-- response-dephasing:start -->
| Gas | Condition | Measurement | Measured (ps) | Response (ps) | Response / measured - 1 |
|---|---|---|---|---|---|
| O2 | S branch, 5.1 atm, revivals to 1.25 periods | chen_2007/T2_O2 | 23.2 | 23.5 | +1% |
<!-- response-dephasing:end -->

**Fused silica.** Stolen & Ippen's gain of fused quartz at 526 nm and 330 cm⁻¹,
measured in a fibre amplifier and from the spontaneous cross section (TC-693);
Hellwarth, Cherlow & Yang's nuclear share of the long-pulse $n_2$ from spontaneous
Raman spectra (TC-694); and the shape of the gain below 25 cm⁻¹, Stolen et al.'s
Eq. (10):

<!-- silica-validation:start -->
| Model | Quantity | Condition | Measurement | Measured | Model | Model / measured - 1 |
|---|---|---|---|---|---|---|
| hollenbeck_2002 | gain | 526.0 nm, 330 cm-1 | stolen_1973/gain_amplifier | 1.5e-11 +- 1.5e-12 | 1.35e-11 | -10% |
| lin_2006 | gain | 526.0 nm, 330 cm-1 | stolen_1973/gain_amplifier | 1.5e-11 +- 1.5e-12 | 1.38e-11 | -8% |
| blow_1989 | gain | 526.0 nm, 330 cm-1 | stolen_1973/gain_amplifier | 1.5e-11 +- 1.5e-12 | 1.46e-11 | -2% |
| hollenbeck_2002 | Raman share of n2 | fused quartz, linear polarisation | hellwarth_1975/nuclear_fraction_SiO2 | 0.21 +- 0.03 | 0.219 | +4% |
| lin_2006 | Raman share of n2 | fused quartz, linear polarisation | hellwarth_1975/nuclear_fraction_SiO2 | 0.21 +- 0.03 | 0.245 | +17% |
| blow_1989 | Raman share of n2 | fused quartz, linear polarisation | hellwarth_1975/nuclear_fraction_SiO2 | 0.21 +- 0.03 | 0.222 | +6% |
| hollenbeck_2002 | gain / peak gain | 5 cm-1 | stolen_1989/low_frequency_gain | 0.004 | 0.0097 | +142% |
| hollenbeck_2002 | gain / peak gain | 10 cm-1 | stolen_1989/low_frequency_gain | 0.009 | 0.0197 | +119% |
| hollenbeck_2002 | gain / peak gain | 20 cm-1 | stolen_1989/low_frequency_gain | 0.026 | 0.0415 | +60% |
| lin_2006 | gain / peak gain | 5 cm-1 | stolen_1989/low_frequency_gain | 0.004 | 0.00531 | +33% |
| lin_2006 | gain / peak gain | 10 cm-1 | stolen_1989/low_frequency_gain | 0.009 | 0.0121 | +34% |
| lin_2006 | gain / peak gain | 20 cm-1 | stolen_1989/low_frequency_gain | 0.026 | 0.0332 | +28% |
| blow_1989 | gain / peak gain | 5 cm-1 | stolen_1989/low_frequency_gain | 0.004 | 0.00528 | +32% |
| blow_1989 | gain / peak gain | 10 cm-1 | stolen_1989/low_frequency_gain | 0.009 | 0.0106 | +17% |
| blow_1989 | gain / peak gain | 20 cm-1 | stolen_1989/low_frequency_gain | 0.026 | 0.0212 | -19% |
<!-- silica-validation:end -->

- **The size.** The Kerr fit takes Stolen & Ippen's gain as data (TC-696), so the
  default's $f_R$ (0.219 ± 0.026) is set by it: the gain comes out 10 % (one standard
  uncertainty) below the measurement, where the $n_2$ data pull, and the share 4 % above
  Hellwarth et al.'s independent one. Its peak at 526 nm, 1.91e-11 cm/W, is 3 % above
  the 1.86e-11 that Stolen & Ippen's cross section implies. Before (the prior 0.18 ±
  0.03, which was derived from these same measurements), $f_R$ was 0.151 and the gain
  40 % low (TC-695). Lin & Agrawal's and Blow & Wood's forms, whose $f_R$ were chosen to
  reproduce gains, are 2-8 % below.
- **The shape below 25 cm⁻¹.** Hollenbeck & Cantrell's gain there is 1.6-2.4 times Stolen
  et al.'s law, which the paper states for their measured spectrum; Lin & Agrawal's
  Boson peak comes closer (1.3 times).


## Response open questions

- **Silica's gain: one measurement.** The Kerr fit's $f_R$ rests on Stolen & Ippen's
  single gain, of a fibre with a fused-quartz core (TC-696). Later fibre gains
  (Dougherty et al. 1995, 1.2e-11 cm/W at 795.5 nm, which Lin & Agrawal quote) are of
  telecommunication fibres, not curated.
- **Silica's gain below 25 cm⁻¹.** Hollenbeck & Cantrell's shape is 1.6-2.4 times Stolen
  et al.'s low-frequency law there, Lin & Agrawal's 1.3 times: it matters for the
  frequency shift of picosecond solitons.
- **Line mixing** of the Q branches of N2, O2 and CO2 above 2-4 amagat (TC-527).
- **The vibrational O and S branches,** 3-19 % of the Q branches.
- **The anisotropic Q lines' widths.** They are taken as the isotropic lines', although
  reorientation also relaxes them.
- **The low densities** below the Dicke cutoff, where the lines are not Lorentzian.
- The gas-side items of the line validation above: the H2 and D2 gain excess, the CO2 and N2O
  stopgaps, N2O's rotational constants.
