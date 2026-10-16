# Third-order support and resonance limits

The model’s data support and the resonance structure are separate constraints.
Read {doc}`nonlinear-conventions` and {doc}`nonlinear-models` first. This chapter
describes the recorded facts; {doc}`nonlinear-extrapolation` examines choices beyond
the fitted data. The {doc}`/nonlinear` examples show the evaluation interface.

## Where the fits are supported

A propagation code that evaluates gamma for every four-wave quartet
$(-\nu_\sigma; \nu_1, \nu_2, \nu_3)$, as the [time-domain form](/user-guide/nonlinear-models.md#the-time-domain-form) does through filters
on each field, meets quartets that no fitted value resembles: a deep-ultraviolet
four-wave-mixing product, or the self-phase modulation of a vacuum-ultraviolet
dispersive wave. The declared range (a fundamental wavelength) and the largest fitted
$W_2$ do not say which of them the fit supports. Three things do: the legs of the data,
the resonances of the gas, and how the polynomial behaves between the two. The tables
of this section are computed by `tools/fit/nonlinear_support.py` from the stored fits
and the curated datasets (TC-864 to TC-868); the page {doc}`/nonlinear_accuracy`
evaluates bounds on the legs against them.

### The legs of the data

Write each quartet as four signed wavenumbers that sum to zero, the output with a minus
sign: the Kerr effect is $(-\nu; \nu, \nu, -\nu)$, the third harmonic
$(-3\nu; \nu, \nu, \nu)$, ESHG $(-2\nu; \nu, \nu, 0)$, the dc Kerr effect
$(-\nu; \nu, 0, 0)$, CARS $(-(2\nu - \nu_s); \nu, \nu, -\nu_s)$ and the pump-probe index
$(-\nu_r; \nu_r, \nu, -\nu)$. $W_2$ is the sum of their squares. For each gas the
largest single leg among the fitted values, the value that sets it, and how many fitted
values have their largest leg (or their $W_2$) within 10 % of the largest:

<!-- nonlinear-support:start -->
| Gas | Largest fitted leg, nm (cm⁻¹) | Set by | W2 there (1e9 cm⁻²) | W4 there (1e18 cm⁻⁴) | Near the edge: leg / W2 / values | Largest fitted W2 (1e9 cm⁻²) | Its Kerr diagonal (nm) | Largest validator leg (nm) |
|---|---|---|---|---|---|---|---|---|
| He | 151.9 (65842) | `bishop_1989/omega_thg_zzzz`, THG at 455.6 nm | 5.780 | 19.490 | 1 / 1 / 52 | 5.780 | 263.1 | 257.3, `shelton_1982/text_He_absolute` (ESHG) |
| Ne | 229.0 (43678) | `hattig_1998/ccsd_eshg_Ne`, ESHG at 457.9 nm | 2.862 | 4.094 | 4 / 1 / 26 | 2.862 | 373.9 | 231.4, `ward_1969/chi_i_Ne` (THG) |
| Ar | 229.0 (43662) | `shelton_1990/table3_Ar_literature`, ESHG at 458.1 nm | 2.860 | 4.089 | 6 / 2 / 41 | 2.860 | 374.0 | 231.4, `ward_1969/chi_i_Ar` (THG) |
| Kr | 244.1 (40974) | `shelton_1990/table3_Kr_literature`, ESHG at 488.1 nm | 2.518 | 3.171 | 2 / 1 / 6 | 2.518 | 398.5 | 231.4, `ward_1969/chi_i_Kr` (THG) |
| Xe | 244.1 (40974) | `shelton_1990/table3_Xe_literature`, ESHG at 488.1 nm | 2.518 | 3.171 | 2 / 1 / 10 | 2.518 | 398.5 | 231.4, `ward_1969/chi_i_Xe` (THG) |
| H2 | 229.0 (43662) | `shelton_1990/table3_H2_literature`, ESHG at 458.1 nm | 2.860 | 4.089 | 6 / 2 / 56 | 2.860 | 374.0 | 231.4, `ward_1969/chi_i_H2` (THG) |
| D2 | 231.5 (43206) | `bishop_1991/thg_par_elec_D2`, THG at 694.4 nm | 2.489 | 3.614 | 5 / 5 / 46 | 2.641 | 389.2 | 435.7, `lundeen_1983/fwm_532_D2` (CARS) |
| N2 | 229.0 (43678) | `hattig_1998/ccsd_eshg_N2`, ESHG at 457.9 nm | 2.862 | 4.094 | 7 / 3 / 37 | 2.862 | 373.9 | 231.4, `ward_1969/chi_i_N2` (THG) |
| O2 | 238.3 (41962) | `shelton_1990/table3_O2_literature`, ESHG at 476.6 nm | 2.641 | 3.488 | 4 / 3 / 20 | 2.641 | 389.2 | 347.1, `ward_1979/chi3_O2` (ESHG) |
| CO2 | 231.5 (43193) | `sekino_1993/thg_CO2`, THG at 694.6 nm | 2.487 | 3.609 | 2 / 2 / 11 | 2.518 | 398.5 | 231.4, `ward_1969/chi_i_CO2` (THG) |
| CH4 | 238.3 (41962) | `shelton_1990/table3_CH4_literature`, ESHG at 476.6 nm | 2.641 | 3.488 | 5 / 4 / 13 | 2.641 | 389.2 | 158.8, `hattig_1998/ccsd_thg_CH4` (THG) |
| SF6 | 229.0 (43662) | `shelton_1990/table3_SF6_literature`, ESHG at 458.1 nm | 2.860 | 4.089 | 2 / 1 / 8 | 2.860 | 374.0 | 231.3, `ward_1984/thg_SF6` (THG) |
| CO | 231.5 (43193) | `sekino_1993/thg_CO`, THG at 694.6 nm | 2.487 | 3.609 | 1 / 1 / 6 | 2.487 | 401.0 | 347.1, `ward_1979/chi3_CO` (ESHG) |
| N2O | 580.0 (17241) | `wahlstrand_2012/n2_N2O`, pump-probe at 800.0 nm | 0.907 | 0.226 | 1 / 1 / 1 | 0.907 | 664.1 | none |
| O3 | 694.6 (14398) | `naves_2011/gamma_kerr_ccsd_t_expt`, Kerr at 694.6 nm | 0.829 | 0.172 | 1 / 1 / 4 | 0.829 | 694.6 | none |
<!-- nonlinear-support:end -->

The largest fitted leg is, for every gas with ESHG or THG data, a harmonic leg: Shelton's
visible ESHG (harmonic at 229-244 nm) or a calculated THG at 694 nm (harmonic at
231.5 nm). It lies 1.6 to 1.7 times further out in wavenumber than the Kerr diagonal of
the same largest $W_2$ (374-401 nm): a cap at the Kerr diagonal cuts into the data. Only
helium's calculation reaches further (151.9 nm). N2O's single pump-probe value reaches
580 nm and O3's calculation 694.6 nm. The validators are evidence, not fitted; Hättig &
Jørgensen's CH4 THG reaches 158.8 nm.

The four legs sum to zero, so each is minus the sum of the other three, and by the
Cauchy-Schwarz inequality $\nu_j^2 \le \tfrac{3}{4} W_2$: no quartet within the
largest fitted $W_2$ has a leg beyond $\sqrt{3 W_{2,\max} / 4}$, which is the
third-harmonic output at that $W_2$.

### Resonance limits

In the sum over states for gamma (Shelton & Rice 1994, eq. 14 on p. 8, and eq. 17 for
the static limit), each term has three energy denominators. With the four frequencies permuted
among them, the first and last intermediate states, which are dipole-coupled to the
ground state, are reached by a single leg (the output's denominator is the three-photon
resonance of the inputs), and the middle state, two dipole steps from the ground state,
by the sum of a pair of legs $|\nu_a + \nu_b|$ (two-photon absorption, or an electronic
Raman resonance when the two legs have opposite signs). Hence two limits per gas:

- **$E_1$**, the lowest level dipole-coupled to the ground state: no leg may reach it;
- **$E_2$**, the lowest level two dipole steps away: no pair sum may reach it.

Ionisation (one leg, a pair, or the output as three photons beyond the ionisation
energy) opens a continuum and makes gamma complex, but $E_1$ and $E_2$ come first.

For an atom with a $J = 0$ ground level, $E_1$ is the lowest odd $J = 1$ level and $E_2$
the lowest even level of $J = 0$, 1 or 2. The two-photon matrix element to an even level,
$\sum_m \langle n|\mathbf d|m\rangle\langle m|\mathbf d|g\rangle$ with the photons'
polarizations and the denominators $1/(\omega_{mg} - \omega_a)$ and
$1/(\omega_{mg} - \omega_b)$, reaches $J = 0$ and 2 through its symmetric part and
$J = 1$ only through its antisymmetric part. The antisymmetric part is proportional to
$\mathbf e_a \times \mathbf e_b$ and to the difference of the two denominators: it
vanishes for two photons of the same frequency, and for any two photons of one linear
polarization, as in a scalar field. So $E_2$ is the lowest even $J = 0$ or 2 level for a
scalar field and the lowest even level of any $J \le 2$ for a vector one (the column *any
polarisation*). Helium's levels are LS-coupled and its triplets below the singlets are
spin-forbidden. A molecule without an inversion centre (CO, N2O, CH4, O3) has states
that are both one- and two-photon allowed, so $E_2 \le E_1$ there.

These limits are for the electronic response in the electric-dipole approximation,
from the ground state; Raman resonances, linewidths, thermally populated states and
weak forbidden transitions are not covered, and being below a level does not establish
an accuracy or negligible absorption. Per gas, $E_1$ and $E_2$ in cm⁻¹, with the
wavelength of the single leg at $E_1$ and of the fundamental whose degenerate pair
reaches $E_2$; the ionisation energy; and $\min(E_1, E_2/2)$ (*Quantities derived from
the facts*, below). Where a polyatomic molecule's one-photon limit is only evidence
about absorption, not a level, the table says what the evidence is:

<!-- nonlinear-resonances:start -->
| Gas | E1, cm⁻¹ (one leg, nm) | E2, cm⁻¹ (a degenerate pair, nm) | E2, any polarisation | Ionisation, cm⁻¹ (one leg, three photons, nm) | Band edge min(E1, E2/2), nm: one polarisation / any | Datasets |
|---|---|---|---|---|---|---|
| He | 171135 (58.4) | 166277 (120.3) | 166277 (120.3) | 198311 (50.4, 151.3) | 120.3 / 120.3 | E1: `nist_asd_2026/level_He_1s2p_1P1`; E2: `nist_asd_2026/level_He_1s2s_1S0`; E2_any: `nist_asd_2026/level_He_1s2s_1S0`; IE: `nist_asd_2026/ie_He_I` |
| Ne | 134459 (74.4) | 149824 (133.5) | 148258 (134.9) | 173930 (57.5, 172.5) | 133.5 / 134.9 | E1: `nist_asd_2026/level_Ne_odd_J1`; E2: `nist_asd_2026/level_Ne_even_J0 + nist_asd_2026/level_Ne_even_J2`; E2_any: `nist_asd_2026/level_Ne_even_J0 + nist_asd_2026/level_Ne_even_J1 + nist_asd_2026/level_Ne_even_J2`; IE: `nist_asd_2026/ie_Ne_I` |
| Ar | 93751 (106.7) | 105617 (189.4) | 104102 (192.1) | 127110 (78.7, 236.0) | 189.4 / 192.1 | E1: `nist_asd_2026/level_Ar_odd_J1`; E2: `nist_asd_2026/level_Ar_even_J0 + nist_asd_2026/level_Ar_even_J2`; E2_any: `nist_asd_2026/level_Ar_even_J0 + nist_asd_2026/level_Ar_even_J1 + nist_asd_2026/level_Ar_even_J2`; IE: `nist_asd_2026/ie_Ar_I` |
| Kr | 80917 (123.6) | 92307 (216.7) | 91169 (219.4) | 112914 (88.6, 265.7) | 216.7 / 219.4 | E1: `nist_asd_2026/level_Kr_odd_J1`; E2: `nist_asd_2026/level_Kr_even_J0 + nist_asd_2026/level_Kr_even_J2`; E2_any: `nist_asd_2026/level_Kr_even_J0 + nist_asd_2026/level_Kr_even_J1 + nist_asd_2026/level_Kr_even_J2`; IE: `nist_asd_2026/ie_Kr_I` |
| Xe | 68045 (147.0) | 78120 (256.0) | 77269 (258.8) | 97834 (102.2, 306.6) | 256.0 / 258.8 | E1: `nist_asd_2026/level_Xe_odd_J1`; E2: `nist_asd_2026/level_Xe_even_J0 + nist_asd_2026/level_Xe_even_J2`; E2_any: `nist_asd_2026/level_Xe_even_J0 + nist_asd_2026/level_Xe_even_J1 + nist_asd_2026/level_Xe_even_J2`; IE: `nist_asd_2026/ie_Xe_I` |
| H2 | 90203 (110.9) | 99165 (201.7) | 99165 (201.7) | 124417 (80.4, 241.1) | 201.7 / 201.7 | E1: `nist_webbook_2026/state_H2_B`; E2: `nist_webbook_2026/state_H2_E + nist_webbook_2026/state_H2_B`; E2_any: `nist_webbook_2026/state_H2_E + nist_webbook_2026/state_H2_B`; IE: `liu_2009/ie_para_H2` |
| D2 | 90634 (110.3) | 99462 (201.1) | 99462 (201.1) | 124745 (80.2, 240.5) | 201.1 / 201.1 | E1: `nist_webbook_2026/state_D2_B`; E2: `nist_webbook_2026/state_D2_E + nist_webbook_2026/state_D2_B`; E2_any: `nist_webbook_2026/state_D2_E + nist_webbook_2026/state_D2_B`; IE: `liu_2010/ie_ortho_D2` |
| N2 | 100817 (99.2) | 68952 (290.1) | 68952 (290.1) | 125667 (79.6, 238.7) | 290.1 / 290.1 | E1: `nist_webbook_2026/state_N2_b + nist_webbook_2026/state_N2_a + nist_webbook_2026/state_N2_ap`; E2: `nist_webbook_2026/state_N2_a + nist_webbook_2026/state_N2_ap`; E2_any: `nist_webbook_2026/state_N2_a + nist_webbook_2026/state_N2_ap`; IE: `trickl_1989/ie_N2` |
| O2 | 49358 (202.6) | 65670 (304.6) | 65670 (304.6) | 97348 (102.7, 308.2) | 304.6 / 304.6 | E1: `nist_webbook_2026/state_O2_B`; E2: `nist_webbook_2026/state_O2_C`; E2_any: `nist_webbook_2026/state_O2_C`; IE: `tonkyn_1989/ie_O2` |
| CO2 | not curated | not curated | not curated | 111127 (90.0, 270.0) | n/a / n/a | IE: `wang_1988/adiabatic_CO2#X` |
| CH4 | none detected to 197.7 nm | not curated | not curated | 101707 (98.3, 295.0) | n/a / n/a | E1: `ityaksov_2008/rayleigh_CH4`; IE: `berkowitz_1987/adiabatic_CH4` |
| SF6 | none detected to 197.7 nm | not curated | not curated | 123806 (80.8, 242.3) | n/a / n/a | E1: `ityaksov_2008/rayleigh_SF6`; IE: `frost_1967/first_threshold_SF6` |
| CO | 64748 (154.4) | 64748 (308.9) | 64546 (309.9) | 112998 (88.5, 265.5) | 308.9 / 309.9 | E1: `nist_webbook_2026/state_CO_A`; E2: `nist_webbook_2026/state_CO_A`; E2_any: `nist_webbook_2026/state_CO_A + nist_webbook_2026/state_CO_I`; IE: `potts_1974/hei_CO#3sigma` |
| N2O | detected in a curated window ending at 45974 (217.5) | not curated | not curated | 103933 (96.2, 288.6) | n/a / n/a | E1: `chan_1994/intervals_present`; IE: `berkowitz_1977/adiabatic_N2O` |
| O3 | detected in a curated window ending at 9553 (1046.8) | not curated | not curated | 102110 (97.9, 293.8) | n/a / n/a | E1: `serdyuchenko_2014/this_work_223`; IE: `mahla_2025/expt_willitsch_O3#2B2 (4b2^-1)` |
<!-- nonlinear-resonances:end -->

- **The atoms:** NIST ASD's levels (`nist_asd_2026`, TC-862). For helium these include the two-photon 1s2s ¹S resonance at a 120.3 nm
  fundamental, the three-photon 1s2p ¹P and 1s3p ¹P resonances at 175.3 and 161.1 nm,
  three-photon ionisation below 151.3 nm. The 1s3p level
  (`nist_asd_2026/level_He_1s3p_1P1`) is kept as that check; it is not a limit, 1s2p
  lying below it. Xenon's lowest even levels give two-photon
  resonances at 256.0 nm (6p[5/2]₂) and 258.8 nm (6p[1/2]₁) as well as the 249.6 nm
  6p[1/2]₀ resonance at which Lehmberg et al. measured a negative n2 (TC-854).
- **The diatomic molecules:** Huber & Herzberg's states (`nist_webbook_2026`, TC-863),
  each from its ground vibrational level. N2's a ¹Πg state, two-photon allowed, sets its
  pair limit at a 290 nm fundamental, far below its dipole-allowed b ¹Πu state (99.2 nm).
  O2's C ³Πg sets its pair limit at 304.6 nm; its Schumann-Runge B state (202.6 nm) is the
  one-photon limit. O2's a ¹Δg and b ¹Σg⁺ states (7882 and 13121 cm⁻¹) are reached by
  pair differences inside the range, but from the triplet ground state they are
  spin-forbidden (the 1.27 µm and 762 nm bands are magnetic-dipole transitions) and are
  not counted.
- **The polyatomic molecules:** the curated absorption data give evidence, not levels.
  Ityaksov et al. detected no absorption by CH4 or SF6 at their sampled wavelengths down
  to 197.7 nm, which is not a proven bound on the lowest level. 5.70 eV (217.5 nm) is
  the lower boundary of Chan et al.'s integration over N2O's B ¹Δ band. CO2's weak
  absorption below its first band near 147.5 nm is not curated. No two-photon state of
  these molecules is curated (TC-869).
- **O3 is the resonant case.** It absorbs throughout its declared range: the Chappuis
  and Wulf bands extend to the end of the curated table at 1046.8 nm (Serdyuchenko et
  al. measured 213-1100 nm; 1046.8 nm is not an onset), weakly (Inn &
  Tanaka put the continuum under the Chappuis bands at an oscillator strength of about
  2e-5), and the strong Huggins and Hartley bands lie at 200-337 nm (Mason et al.'s
  intervals), where the third harmonic of its range falls. Its fit is a non-resonant
  calculation (Naves et al.).

### Does gamma increase with each frequency?

**The fits do, everywhere.** The derivative of the model with respect to one leg's
$\nu_j^2$ is $\gamma_0 (A + 2B W_2 + 2B'\nu_j^2)$. Every fit has $A > 0$, $B \ge 0$ and
$B + B' > 0$, and $\nu_j^2 \le W_2$, so the derivative is at least $\gamma_0 A > 0$ at any
$W_2$, within the fitted data and beyond (the minimum, $A$, at the static limit).

**The physical gamma does in every curated case, but no theorem says so.** Shelton &
Rice's eq. 17 shows that the sum over states has a negative term, so that even the static
gamma could in principle be negative (a two-level system), though no such case had been
reported (p. 10). Every curated dispersion curve rises from the static limit toward the
first resonance: their Fig. 1 for the hydrogen atom, Bishop & Pipin's helium, Hättig &
Jørgensen's Ne, Ar, N2 and CH4, Bishop et al.'s H2 and D2. Farther beyond the fitted
data the calculations generally lie above the polynomial, though not at every point near
the window edge (argon's THG at 632.8 nm sits at 1.001; *Third-order accuracy*); this
covers particular processes and gases and proves no general lower bound (TC-869). The
polynomial's increase is an algebraic property on non-negative squared legs, capped
ones included, not a physical continuation of a quartet. The statement
fails near and above a resonance: above a two-photon resonance gamma can turn over and
change sign, as xenon's n2 does at 248.4 nm (TC-854).

### Quantities derived from the facts

Two quantities follow from the facts above by definition:

- **$\sqrt{3 W_{2,\max}/4}$**, the largest leg any quartet within the largest fitted
  $W_2$ can have. A leg filtered to $\min(|\nu_j|, \sqrt{3 W_{2,\max}/4})$ leaves every
  fitted value and every quartet within the fitted $W_2$ as the fit gives them.
- **$\min(E_1, E_2/2)$**, a sufficient bound: if every leg of a quartet, the generated
  output included, is below it, no leg reaches $E_1$ and no pair sum ($\le$ twice the
  largest leg) reaches $E_2$. A bound on a field's input band alone does not bound the
  outputs it generates, which must be included. For one linear polarisation the SPM
  quartet $(-\nu, \nu, \nu, -\nu)$ reaches the scalar $E_2$ at the edge; for any
  polarisation the edge need not be an SPM resonance, since the atoms' even $J = 1$
  level is not reached by two photons of equal frequency. Exceeding it does not prove a
  quartet resonant, particularly in a sparse spectrum; the record gives no linewidths
  or tensor amplitudes. It is not defined where $E_1$ or $E_2$ is not curated, or where
  $E_1$ is only absorption evidence (TC-869).

How a propagation code uses them, and what it does beyond the data, is its own
modelling decision (TC-868). The page
{doc}`/nonlinear_accuracy` gives the evidence such a decision can draw on: per-leg caps
evaluated against the data, the quartets within the fitted $W_2$, a Kerr carrier and
representative four-wave-mixing quartets, and the polynomial against the calculations beyond
the data.

<!-- nonlinear-derived:start -->
| Gas | Largest leg within the fitted W2, √(3 W2,max / 4), nm (cm⁻¹) | min(E1, E2/2), nm: one polarisation / any | Largest fitted pair, as a degenerate fundamental (nm) | Note |
|---|---|---|---|---|
| He | 151.9 (65842) | 120.3 / 120.3 | 455.6 | - |
| Ne | 215.9 (46327) | 133.5 / 134.9 | 457.9 | - |
| Ar | 215.9 (46311) | 189.4 / 192.1 | 458.1 | - |
| Kr | 230.1 (43459) | 216.7 / 219.4 | 488.1 | - |
| Xe | 230.1 (43459) | 256.0 / 258.8 | 488.1 | - |
| H2 | 215.9 (46311) | 201.7 / 201.7 | 458.1 | - |
| D2 | 224.7 (44507) | 201.1 / 201.1 | 476.6 | - |
| N2 | 215.9 (46327) | 290.1 / 290.1 | 457.9 | - |
| O2 | 224.7 (44507) | 304.6 / 304.6 | 476.6 | - |
| CO2 | 230.1 (43459) | n/a | 488.1 | E1 not curated; E2 not curated |
| CH4 | 224.7 (44507) | n/a | 476.6 | E1: absorption evidence only; E2 not curated |
| SF6 | 215.9 (46311) | n/a | 458.1 | E1: absorption evidence only; E2 not curated |
| CO | 231.5 (43193) | 308.9 / 309.9 | 694.3 | - |
| N2O | 383.4 (26082) | n/a | 672.5 | E1: absorption evidence only; E2 not curated |
| O3 | 401.0 (24937) | n/a | 694.6 | E1: absorption evidence only; absorbs inside the declared range; E2 not curated |
<!-- nonlinear-derived:end -->

### The support record

Each model carries these facts as a public record, `support`
({class}`~optmat.gas.nonlinear.HyperpolarizabilitySupport`), so that a propagation code
reads them instead of copying this page's tables
({doc}`OEP 5 </oeps/archive/oep-0005-hyperpolarizability-support>`). Wavenumbers are inverse
vacuum wavelengths in m⁻¹, $W_2$ is in m⁻², and every fact carries its sources:

- `max_w2`, `leg_max` and `pair_max`: the largest fitted $W_2$, single leg and pair
  sum, over the points the fit admits after its windows and row selections; validators
  and windowed-out values do not count. Each names its *witnesses*, the fitted points
  that attain it, with the dataset, row, process, wavelengths and all four signed legs.
  These are extents of the sampled data, not an interpolation domain or an accuracy.
- `e1`, `e2` and `e2_any`: the resonance limits of the table above, for one linear
  polarisation and for any. Each is a
  {class}`~optmat.gas.nonlinear.SupportLimit` with its `kind` (`level`, `band_origin`,
  `absorption_evidence`, ...), its datasets and where in them, how it is derived, its
  register entries, its qualifications and its uncertainty as the source states it, or
  the fact that it states none. Where a limit is the smallest of several candidates, the
  candidates are listed. A limit that is not curated is None, with its reason in
  `missing_reasons`; absorption evidence (CH4, SF6, N2O, O3) records what was measured
  and is never taken as a level.
- `ionisation`: the ionisation energy over $hc$, the dataset of
  {mod}`optmat.gas.ionisation` (ozone's: the lower of Mahla et al.'s thresholds).
- Derived: `leg_max_in_w2`, $\sqrt{3 W_{2,\max}/4}$, and
  `band_edge(polarisation="linear")`, $\min(E_1, E_2/2)$ with `e2` (`"linear"`) or
  `e2_any` (`"any"`): None where either limit is missing or is not a level or a band
  origin.

The record states facts and carries no rule (TC-868). It changes no evaluated gamma and
is not part of the model's digest; it has its own, which a propagation code stores
beside the model's. Its `model_digest` is the digest of the model it describes: a model
built from a caller's parameters has no record (None), and a curated model whose record
is missing or stale raises a {class}`~optmat.exceptions.DataError`. The records are
derived from the stored fits, the curated datasets and the selection
`data/gas/nonlinear_support.toml` (which names each limit's datasets and why), and
stored in `data/gas/models/nonlinear_support/<gas>.json`; a test re-derives them
(TC-870).

```python
from optmat.gas.gases import Gas
from optmat.gas.nonlinear import hyperpolarizability_models

argon = hyperpolarizability_models()[Gas.Ar]
support = argon.support
support.leg_max.wavenumber, support.leg_max.witnesses[0].process  # 4.37e6 m^-1, eshg
support.e1.datasets  # ('nist_asd_2026/level_Ar_odd_J1',)
support.e2_any.candidates  # the even J = 0, 1 and 2 levels; J = 1 is the smallest
support.band_edge(), support.band_edge("any")  # 189.4 and 192.1 nm, in m^-1
support.digest  # its own identity, beside argon.identity.digest
assert hyperpolarizability_models()[Gas.SF6].support.band_edge() is None
```
