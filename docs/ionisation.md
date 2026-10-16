# Ionisation

A strong-field ionisation rate (ADK, PPT, molecular ADK) is evaluated from a few
properties of the atom or molecule. optmat supplies those properties, curated and with
their sources (TC-439, TC-444). For the noble gases it also gives the rates
themselves, validated against calculations and measurements: {doc}`ionisation_rates`
(TC-652). For the molecules it gives a rate set by all-electron calculations, and the
forms that propagation codes use: {doc}`molecular_ionisation_rates` (TC-673, TC-677).
Everything here comes from {func}`optmat.gas.ionisation.ionisation_parameters`:

```python
from optmat.gas.ionisation import ionisation_parameters

argon = ionisation_parameters("Ar")
argon.ionization_energy  # J, NIST ASD
argon.n_star, argon.l_star  # the effective quantum numbers of ADK and PPT
argon.hartree_coefficient_squared  # |C_n*l*|^2, Ammosov et al. eq. (19)
argon.structure["tong_2002"].coefficients  # C_l of the orbital, Tong et al. 2002
argon.polarizability, argon.ion_polarizability  # C m^2/V, for the Stark shift
argon.barrier_suppression_alpha  # Tong & Lin 2005
argon.barrier_suppression_zhang  # Zhang, Lan & Lu 2014: a1, a2, a3
argon.reasons  # why each value was chosen
```

## What is provided

- **The first ionisation energy.** The atoms' come from the NIST Atomic Spectra
  Database (TC-443). The molecules' come from their primary determinations (TC-542):
  - H2 and D2: Liu et al. (2009, 2010), N = 0 → N⁺ = 0;
  - N2: Trickl et al. 1989;
  - O2: Tonkyn et al. 1989;
  - CH4: Berkowitz et al.'s 1987 photoion onset;
  - SF6: Frost et al.'s 1967 first threshold, whose kind the paper does not state;
  - N2O: Berkowitz & Eland 1977;
  - CO2: Wang et al. 1988.

  The vertical energies of N2O and CO2 are Turner et al.'s 1970 He I band maxima as
  Potts & Williams 1974 reprint them (the book itself is not in the library), and CH4's
  is Heilbronner et al.'s 1977 first Jahn-Teller maximum (TC-547). The WebBook
  (`nist_webbook_2026`) is no longer the source of any recommended value. Its eV values
  for H2, D2 and O2 are the cm⁻¹ values converted with the 1973 CODATA factor, 8 ppm
  high (TC-505).
  - The recommended energy is the adiabatic one.
  - Where the vertical energy differs it is given too (O2, CH4, SF6, N2O, CO2). A
    tunnelling rate at fixed nuclei sees the vertical energy, and the difference, 1 eV
    for CH4, enters the rate's exponent.
- **The ionising orbital:** its name, an atom's valence l, and the magnetic quantum
  number m along a molecule's axis. For an atom's p shell, m = 0 is the member that
  ionises fastest in a linearly polarised field (Ammosov, Delone & Krainov 1986,
  eq. 14).
- **The structure coefficients** C_l of the orbital's asymptotic tail, in the
  normalisation of Tong, Zhao & Lin (2002), eqs. (1)-(2):
  Ψ = Σ_l C_l F_l(r) Y_lm, with F_l → r^(Z/κ − 1) e^(−κr) and κ = (2 I_p)^(1/2) in
  atomic units.
  - Every curated tabulation is available by name, with a default. The molecules
    default to Zhao et al. 2010 and the atoms to Tong et al. 2002.
  - CH4 and SF6 have Zhao et al.'s 2011 coefficients, combined across degenerate
    orbitals; N2O has Zhao et al.'s 2017 table (TC-681, TC-684, TC-675). These are
    available by name, but `default_structure` is unset because their independent
    validation is limited. The molecular rate models state which set they use.
- **The derived quantities of ADK and PPT:** n* = Z/κ (Ammosov et al., eq. 2) and
  l* = n* − 1 (their eq. 17), and the Hartree coefficient |C_n*l*|^2 of their eq. (19).
  {func}`~optmat.gas.ionisation.adk_coefficient_squared` is the Stirling form, their
  eq. (20), which is misprinted in the paper and corrected here (TC-441).
- **The static polarizabilities:**
  - **The neutral's:** the electronic static limit of the gas's fitted dispersion model
    (TC-441), with its uncertainty. The nuclei do not follow an optical field, so the
    infrared (vibrational) part is left out.
  - **The cation's:** Wang et al. 2021's CCSD(T) values for He+, Ne+ and Ar+. The
    others are requested.

  The difference α₀ − α₊ sets the Stark shift (α₀ − α₊)F²/2 of the ionisation energy
  in a static field F.
- **The model constants.** These are fitted inside one rate formula, not properties
  of the gas, and each is labelled as such:
  - Tong & Lin's (2005) barrier-suppression α, fitted for He, Ne and Ar.
  - Zhang, Lan & Lu's (2014) a1, a2, a3 for He to Xe. The rate is
    W_ADK exp(a1 x² + a2 x + a3), with x = F/F_b. Their eq. (8) prints the exponent's
    sign as minus, but their own Table IX needs plus (TC-509, TC-542).
  - Talebpour, Yang & Chin's (1999) effective charges for N2 (0.9) and O2 (0.53).
    These were fitted within PPT with l = m = 0. O2's was fitted with I_p = 12.55 eV
    (O2⁺ v⁺ = 1), not with the adiabatic 12.07 eV (TC-508).

The parameters of every gas:

<!-- ionisation-parameters:start -->
| Gas | I_p (eV) | Source | Vertical (eV) | Orbital | n* | Structure coefficients | alpha (au) | alpha+ (au) | alpha_BSI | Zhang a1, a2, a3 | Z_eff (PPT, l = 0) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| He | 24.5874 | nist_asd_2026/ie_He_I |  | 1s (m = 0) | 0.7439 | tong_2002: C0 = 3.13 | 1.384 | 0.2811 | 7 | 0.1355, -0.8621, 0.021562 |  |
| Ne | 21.5645 | nist_asd_2026/ie_Ne_I |  | 2p (m = 0) | 0.7943 | tong_2002: C1 = 2.1 | 2.661 | 1.311 | 9 | 0.10061, -1.04832, -0.07542 |  |
| Ar | 15.7596 | nist_asd_2026/ie_Ar_I |  | 3p (m = 0) | 0.9292 | tong_2002: C1 = 2.44 | 11.08 | 6.834 | 9 | 0.16178, -1.50441, 0.32127 |  |
| Kr | 13.9996 | nist_asd_2026/ie_Kr_I |  | 4p (m = 0) | 0.9858 | tong_2002: C1 = 2.49 | 16.76 |  |  | 0.1464, -1.36533, 0.02055 |  |
| Xe | 12.1298 | nist_asd_2026/ie_Xe_I |  | 5p (m = 0) | 1.0591 | tong_2002: C1 = 2.57 | 27.08 |  |  | 0.2108, -1.88482, 0.574281 |  |
| H2 | 15.4258 | liu_2009/ie_para_H2 |  | HOMO 1 sigma_g (m = 0) | 0.9392 | zhao_2010: C0 = 1.78, C2 = 0.11, C4 = 0 | 5.433 |  |  |  |  |
| D2 | 15.4665 | liu_2010/ie_ortho_D2 |  | HOMO 1 sigma_g (m = 0) | 0.9379 | zhao_2010: C0 = 1.78, C2 = 0.11, C4 = 0 | 5.364 |  |  |  | 0.72 |
| N2 | 15.5807 | trickl_1989/ie_N2 |  | HOMO 3 sigma_g (m = 0) | 0.9345 | zhao_2010: C0 = 2.68, C2 = 1.1, C4 = 0.06 | 11.74 |  |  |  | 0.9 |
| O2 | 12.0696 | tonkyn_1989/ie_O2 | 12.3070 | HOMO 1 pi_g (m = 1) | 1.0617 | zhao_2010: C2 = 0.52, C4 = 0.03 | 10.59 |  |  |  | 0.53 |
| CH4 | 12.6100 | berkowitz_1987/adiabatic_CH4 | 13.6000 | HOMO 1 t2 (m = 0) | 1.0387 |  | 17.26 |  |  |  |  |
| SF6 | 15.3500 | frost_1967/first_threshold_SF6 | 15.6900 | HOMO 1 t1g (m = 0) | 0.9415 |  | 30.06 |  |  |  |  |
| N2O | 12.8860 | berkowitz_1977/adiabatic_N2O | 12.8900 | HOMO 2 pi (m = 1) | 1.0275 |  | 19.7 |  |  |  |  |
| CO2 | 13.7780 | wang_1988/adiabatic_CO2#X | 13.7900 | HOMO 1 pi_g (m = 1) | 0.9937 | zhao_2010: C2 = 1.97, C4 = 0.4, C6 = 0.04 | 17.53 |  |  |  |  |
<!-- ionisation-parameters:end -->

## How well the parameters give the rates

The parameters are validated against their own primary sources (spectroscopy for the
ionisation energies, quantum chemistry for the coefficients). A rate tests the
parameters and the formula together, so the comparison below is reported, never fitted
to. It is of the **uncorrected** tunnelling rate, to leading order in the field; the
rates of {doc}`ionisation_rates` carry the correction at finite field, and are
compared with the same references there.

The references are the ab initio static-field rates curated so far (tier 1 of TC-439):
- Parker et al.'s full two-electron TDSE for He;
- Majety & Scrinzi's haCC for He, H2 and N2 (tRecX);
- Saenz's frozen-nuclei H2;
- Zhang, Lan & Lu's single-active-electron rates for Ar, in the GSZ model potential.
  This is a one-electron calculation, so it tests the formula, not the
  multi-electron physics.

The table gives the static rate of Tong et al.'s eq. (7) for **one electron** in the
orbital, aligned with the field, divided by the reference, at the fields (au) each
reference gives:

<!-- ionisation-validation:start -->
| Gas | Reference | F = 0.03 | F = 0.045 | F = 0.05 | F = 0.053 | F = 0.06 | F = 0.08 | F = 0.094 | F = 0.1 | F = 0.12 | F = 0.14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| He | Parker et al. 2009, full TDSE |  | 1.03 | 1.01 |  | 0.87 | 0.82 |  | 0.89 | 0.96 | 1.05 |
| He | Majety & Scrinzi 2015, haCC(5) |  |  |  |  |  | 0.83 |  | 0.90 | 0.98 | 1.07 |
| H2 | Saenz 2002, frozen nuclei, parallel |  |  |  | 1.67 | 1.42 | 1.51 |  | 1.87 | 2.22 |  |
| H2 | Majety & Scrinzi 2015, parallel |  |  |  |  | 2.22 | 2.06 |  | 2.21 | 2.49 |  |
| N2 | Majety & Scrinzi 2015, parallel | 1.17 |  | 1.95 |  | 2.66 | 4.99 |  | 8.57 | 12.99 | 17.63 |
| Ar | Zhang et al. 2014, SAE (GSZ potential) |  |  |  | 2.18 |  |  | 3.75 |  |  |  |
<!-- ionisation-validation:end -->

- **Helium:** the one-electron rate with Tong et al.'s C_0 = 3.13 matches the exact
  two-electron rate to 1-3 % at the lowest fields. Between 0.06 and 0.12 au it is 4-18 %
  low, and it rises above the exact rate as the field approaches barrier suppression.
  - Tong & Lin's barrier-suppression factor lowers it everywhere, to 0.63-0.89 of the
    exact rate.
  - So an additional occupancy factor of 2 doubles an
    already correct low-field rate. With the correction at finite field, which the
    rate needs at every field (TC-655), the exact rate is 1.3 to 1.6 times the
    one-electron rate between 0.06 and 0.25 au: neither 1 nor 2 (TC-658).
- **H2:** the one-electron rate is 1.4-2.2 times the frozen-nuclei rate and 2.1-2.5
  times the fixed-nuclei one.
- **N2:** the one-electron rate is 1.2 times the reference at 0.03 au, and the excess
  grows steeply towards N2's barrier-suppression field, κ⁴/16 = 0.08 au.
- **Ar:** 2.2 times the SAE rate at 0.053 au (0.64 of the barrier-suppression field)
  and 3.8 times at 0.094 au. With the same C_l = 2.44, Zhang et al.'s fitted correction
  reproduces their Table IX to within 5 % over 0.64-3.8 of that field (TC-509).

No ab initio multi-electron or calibrated measured rates are curated yet for Ne, Ar,
Kr, Xe, O2, CO2, CH4, SF6 and N2O. Majety & Scrinzi's Ne and Ar are figures only. The
measured yields now in the library are figures only (TC-508, TC-510 to TC-512):
- Larochelle et al. 1998, Guo et al. 1998 and Talebpour et al. 1999 calibrate their
  intensities against a model rate (ADK or PPT), so they are not independent rate
  references.
- Walker et al. 1994 calibrate theirs on xenon photoelectron spectra, to about 30 %.

Where the literature leaves a gap, TDSE reference calculations (tRecX) will fill it, as
validators ranked below published ab initio results.

<span id="against-solprop"></span>

## Source interpretations

The parameters must be read with their definitions. Berkowitz et al.'s rotational
model puts methane's rotationless threshold at about 12.615 eV; SF6 has two
photoelectron onsets. Wang et al.'s Hartree–Fock polarizabilities are 3–11 % below
the correlated values. These distinctions matter when comparing parameter sets.

Talebpour, Yang and Chin fitted O2's effective charge of 0.53 with $I_p = 12.55$ eV.
Pairing it with 12.07 eV gives a faster rate than the fitted one; their $l = 0$ is a
model choice, not O2's $\pi_g$, $m = 1$ orbital. See
{doc}`molecular_ionisation_rates` for the resulting rate models. Zhang et al.'s
coefficients are read from Table XI with the plus sign in eq. (8), which reproduces
their own Table IX (TC-542).

## Open questions

- **Which energy for tunnelling:** adiabatic or vertical. The two differ by 1 eV for
  CH4 and 0.34 eV for SF6 (TC-444).
- **Occupancy:** the factor for doubly occupied or degenerate orbitals that belongs
  with each set of coefficients. The noble-gas defaults use an effective
  electron count with gas-specific evidence and limitations (TC-658, TC-848, TC-850;
  {doc}`user-guide/ionisation-static`). The molecular coefficients and occupancy must
  be considered together ({doc}`molecular_ionisation_rates`).
- **Missing parameters:** polarizabilities of Kr+, Xe+ and the molecular cations.
  Structure coefficients are curated for CH4, SF6 and N2O, but need independent
  validation. The vertical energies of
  N2O and CO2 rest on a reprint of Turner et al. 1970.
- **Rate references:** the current atomic references include static and laser-field
  all-electron calculations ({doc}`user-guide/ionisation-static`,
  {doc}`user-guide/ionisation-laser`). Molecular gaps and the available H2 and N2
  checks are described in {doc}`user-guide/molecular-ionisation-validation`.
