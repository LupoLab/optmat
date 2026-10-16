# Noble-gas ionisation: static references

The static references constrain orbital coefficients, finite-field corrections, and the effective number of electrons. Start with {doc}`ionisation-theory`; {doc}`ionisation-laser` then tests how the atom responds in a laser field. Dataset identifiers below lead to records in the {doc}`/the-data/sources` catalogue.

## The effective number of electrons

In the weak-field limit the rate of a closed shell is that of one electron times the
number of electrons in the orbital, 2 for $m = 0$ (Tolstikhin, Madsen & Morishita 2014,
section II G: "the SAEA result for the ionization rate must be multiplied by 2, due to
the summation over $M_S'$"). At the fields where atoms ionise in a laser pulse the
first-order corrections are not small, and those of the atom exceed those of the
one-electron model. For helium everything is known (Trinh et al. 2015, Table 1): the
exact static rate is 0.74 of the leading order at 0.08 au and 0.43 at 0.20 au.

<!-- ionisation-rate-electrons:start -->
| Gas | N_eff | C_l | alpha | I_p (eV) | I_p to 2P1/2 (eV) | F_b (au) | I at F_b (W/cm2) |
|---|---|---|---|---|---|---|---|
| He | 1.47 | 3.13 | 7 | 24.5874 |  | 0.2041 | 1.46e+15 |
| Ne | 1 | 2.1 | 9 | 21.5645 | 21.6613 | 0.1570 | 8.65e+14 |
| Ar | 1 | 2.44 | 9 | 15.7596 | 15.9371 | 0.0839 | 2.47e+14 |
| Kr | 1.35 | 2.49 | 9 | 13.9996 | 14.6654 | 0.0662 | 1.54e+14 |
| Xe | 1.5 | 2.57 | 9 | 12.1298 | 13.4363 | 0.0497 | 8.66e+13 |
<!-- ionisation-rate-electrons:end -->

- **Helium:** the ratio of Parker et al.'s two-electron static rate to the one-electron
  rate of this model; $N_\text{eff}$ is its geometric mean from 0.07 to 0.25 au. Four
  other all-electron calculations agree with Parker's within 4 % above 0.09 au.
- **Neon and argon:** Jagau's (2016) complex-scaled coupled-cluster rates, with all
  electrons correlated, are tabulated in his supplement and curated (`jagau_2016`,
  TC-668). Over the one-electron rate of this model they are 1.03 to 1.10 for neon
  from 0.16 to 0.24 au, and 0.83 to 1.00 for argon from 0.08 to 0.16 au (TC-669).
  Majety & Scrinzi's (2015) rates, which are in a figure only and are not curated,
  agree (TC-658). For argon the measured yields of Wallace et al.
  ({doc}`ionisation-measured`) support the choice independently.
- **$N_\text{eff}$ is a calibration of this formula, not a count of active
  electrons.** Against a one-electron rate with the exact asymptotic coefficient, the
  many-electron factor of a doubly occupied orbital is 2 in the weak-field limit
  (Tolstikhin, Madsen & Morishita 2014, section II G); the values here, about 1 to
  1.5, are against Tong, Zhao & Lin's coefficients, Tong & Lin's static correction and
  this Coulomb factor, and part of their distance from 2 is the formula's (TC-851).
  What they rule out is a larger count: the occupancy of 2 per $m$ (6 for a $p$ shell)
  that some propagation codes use.
- **The weak-field limit on this formula's scale.** Tolstikhina, Morishita &
  Tolstikhin (2014, Table I) give the atom's coefficient $g = \sqrt{p} R_N g_N$ of
  the many-electron weak-field asymptotic theory from Hartree-Fock orbitals, the
  occupation factor $p = 2$ and the relaxation $R_N$ included. By their Eq. (11) one
  electron of this model has $g_1 = \sqrt{2(2l+1)/\kappa}\, 2^{-1/\kappa} C_l$,
  with which their rate $|g|^2 W(F)$ is this model's static rate at weak field
  (hydrogen: $C_l = 2$, $g = \sqrt 2$, as they print). So $(g/g_1)^2$ is the atom's
  weak-field rate over the model's one electron
  ({func}`~optmat.gas.ionisation_validation.weak_field_factors`; both at the
  Hartree-Fock $\kappa$, 0.8 to 3.6 % above the measured ionisation energies'):

<!-- ionisation-rate-weak-field:start -->
| Gas | kappa (HF, au) | g, the atom (HF) | g, the model's one electron | atom / one electron at weak field | N_eff (default) |
|---|---|---|---|---|---|
| He | 1.3550 | 2.942 | 2.280 | 1.67 | 1.47 |
| Ne | 1.3042 | 1.919 | 2.647 | 0.53 | 1 |
| Ar | 1.0872 | 3.977 | 3.030 | 1.72 | 1 |
| Kr | 1.0239 | 4.091 | 3.063 | 1.78 | 1.35 |
| Xe | 0.9563 | 4.212 | 3.118 | 1.82 | 1.5 |
<!-- ionisation-rate-weak-field:end -->

  For He, Ar, Kr and Xe it is 1.7 to 1.8: the factor 2 of the doubly occupied orbital,
  less the 9 to 17 % by which Tong, Zhao & Lin's $C_l^2$ exceed Hartree-Fock's. At
  working fields the all-electron calculations and Wallace et al.'s yields give less
  (helium 1.47; argon about 1, rising towards 1.7 below 0.07 au, as the next bullet
  shows; krypton's 1.35 and xenon's 1.5 multiply the fine-structure weighting, 0.75 to
  0.85, so about 1.1 without it): the many-electron factor falls with the field from
  its weak-field value, which a constant $N_\text{eff}$ does not follow. Neon's
  Hartree-Fock coefficient is half of what Tong, Zhao & Lin's $C_l$ implies (0.53),
  where Jagau's coupled-cluster rates agree with the formula and 1 electron (1.03 to
  1.10): one of the two coefficients is wrong for neon, and which is not settled
  (TC-852). For helium the accurate two-electron $g$ (2.932) is the Hartree-Fock one
  within 0.4 %.
- **Towards low field the number rises**, as the weak-field limit of 2 requires: for
  argon 1.2 at 0.07 au and 1.6 at 0.06 au (1.7e14 and 1.3e14 W/cm²), for neon 1.2 at
  0.14 au. A constant of 1 is then low by that much. The references are least certain
  there: their rates are below 1e-4 au, neon's are not smooth in the field below
  0.14 au, and argon's at 0.05 au is 3.7 times Majety & Scrinzi's.
- **Argon in a laser pulse** (haCC's helium is 0.76 to 0.81 of the exact two-electron
  yield at 390 nm, TC-1050: the haCC numbers below may be low by as much): haCC's
  all-electron calculation of argon (A2's 6 fs pulse
  at 800 nm; the Scrinzi group's `Ar_qz` data, four or three ion channels), against
  one $m = 0$ electron at the same ionisation energy and with the same measure, the
  norm-loss after the pulse ({func}`~optmat.gas.ionisation_validation.hacc_ratios`):
  0.76, 0.80 and 0.87 at 4.0e13, 9.6e13 and 2.3e14 W/cm² (Keldysh parameter 1.8 to
  0.8), 0.81 with three channels. With all three 3p channels and the $m = \pm 1$
  electrons (mMax 1) the run at 9.6e13 W/cm² gives 0.86 of one electron, about 0.79
  when its angular momentum is converged (TC-847, TC-861). In a 20 fs pulse the
  three-channel set-up gives 0.84 against its 0.81 at 6 fs, against one electron run at
  exactly its ionisation energy in both pulses; with lMax 30, 0.81 against 0.75: the
  ratio rises by 3 to 8 % with the pulse's length, and the angular correction (6.8 %
  at 6 fs, 3.1 % at 20 fs) does not cancel between them (TC-861). With the
  frequency-dependent Coulomb factor, van der Hart's all-electron R-matrix rates at 390 nm give
  about 0.87 (TC-842): in a laser pulse argon ionises like about 0.8 to 0.9 of one
  electron, the default's 1 within 20 %.
- **Published many-electron calculations in short pulses** (TC-858 to TC-860): the
  time-dependent CASSCF of argon at 1200 nm (Pathak et al. 2020) and 800 nm (Sato et
  al. 2018), of krypton at 800 nm (Pathak et al. 2022) and the converged two-electron
  helium at 800 nm (Sato et al. 2016), against one m = 0 electron (Tong and Lin's or
  part G's potential) run in the same sin^2 pulses and read with their measure, the
  probability outside 20 au at the end of the pulse
  ({func}`~optmat.gas.ionisation_validation.published_atom_ratios`; the rows
  `pathak_2020/fig2_sae_Ar` are their own single-electron curve against one electron in
  their potential):

<!-- ionisation-rate-atoms:start -->
| Published | Gas | nm | I (W/cm2) | Published, outside 20 au | One electron (tRecX) | Ratio |
|---|---|---|---|---|---|---|
| `pathak_2020/fig2_tdcasscf_Ar` | Ar | 1200 | 1.00e+14 | 0.0007255 | 0.0006582 | 1.10 |
| `pathak_2020/fig2_tdcasscf_Ar` | Ar | 1200 | 2.00e+14 | 0.03102 | 0.03275 | 0.95 |
| `pathak_2020/fig2_tdcasscf_Ar` | Ar | 1200 | 4.00e+14 | 0.3401 | 0.3917 | 0.87 |
| `pathak_2020/fig2_sae_Ar` | Ar | 1200 | 1.00e+14 | 0.0004512 | 0.000506 | 0.89 |
| `pathak_2020/fig2_sae_Ar` | Ar | 1200 | 2.00e+14 | 0.02225 | 0.02559 | 0.87 |
| `pathak_2020/fig2_sae_Ar` | Ar | 1200 | 4.00e+14 | 0.2918 | 0.33 | 0.88 |
| `pathak_2022/fig2_tdcasscf_Kr` | Kr | 800 | 2.00e+14 | 0.1205 | 0.1331 | 0.91 |
| `sato_2018/fig3a_casscf_Ar_800nm` | Ar | 800 | 6.00e+14 | 0.5078 | 0.6373 | 0.80 |
| `sato_2016/fig6_casscf_He_800nm` | He | 800 | 4.00e+14 | 0.0003029 | 0.0001775 | 1.71 |
| `sato_2016/fig6_casscf_He_800nm` | He | 800 | 8.00e+14 | 0.01589 | 0.009906 | 1.60 |
| `sato_2016/fig6_casscf_He_800nm` | He | 800 | 1.40e+15 | 0.1608 | 0.105 | 1.53 |
<!-- ionisation-rate-atoms:end -->

  Argon is 0.87 to 1.10 of one electron at 1200 nm, falling with intensity, and 0.80 at
  800 nm and 6e14 W/cm2, near saturation, where a ratio of probabilities understates
  one of rates; krypton 0.91 at 2e14 W/cm2; helium 1.53 to 1.71, as Parker et al.'s
  (1.47) and haCC's (1.65). These are ratios to one electron of Tong and Lin's
  potential: Pathak et al.'s own potential (Muller and Kooiman's) gives 0.77 to 0.84 of
  its yield at an ionisation energy 0.3 % higher, so 'one electron' is uncertain by
  that much; their single-electron curve is 0.87 to 0.89 of tRecX's in that potential,
  unexplained. The two krypton and argon calculations at 800 nm assume the group's
  sin^2 envelope on the field, which the papers do not state.
- **Krypton and xenon:** no all-electron rate of either exists. Santra et al.'s (2006)
  krypton rate is an independent-electron sum over the six 4p spin-orbitals in a
  Herman-Skillman potential with the spin-orbit coupling, not an all-electron
  calculation; the default is 1.04 times it at 0.10 au, a check of the model in another
  one-electron potential, not of the number of electrons (TC-851). Krypton's 1.35 and
  xenon's 1.5 rest on Wallace et al.'s yields through the exact one-electron
  calculations in their pulse (TC-846, TC-849, TC-850), and put the three gases at the
  same intensity scale (1.02 to 1.03). For xenon one many-electron calculation in a
  stated pulse is curated: Hartmann et al.'s (2022) R-matrix (RMT) ground-state
  survival in a 5.3 fs pulse at 759 nm, 1.1 to 2.5e14 W/cm² (`hartmann_2022`, TC-1059).
  It counts excitation as well as ionisation; a one-electron calculation in the same
  field, read the same way, is to be compared with it.

## Static rates

Every number is the model's value over the reference's.

<!-- ionisation-rate-static:start -->
| Gas | Reference | Model | 0.05 au | 0.06 au | 0.07 au | 0.08 au | 0.1 au | 0.12 au | 0.15 au | 0.2 au | 0.25 au | 0.3 au |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| He | Parker et al. 2009, two electrons | `default` | 1.25 | 1.05 | 0.92 | 0.93 | 0.92 | 0.97 | 0.99 | 1.10 | 1.16 | 1.17 |
| He | Scrinzi et al. 1999, two electrons | `default` |  |  |  | 1.25 | 0.97 | 0.98 | 1.00 | 1.09 | 1.16 | 1.19 |
| He | Majety & Scrinzi 2015, haCC(5) | `default` |  |  | 0.94 | 0.95 | 0.97 | 0.98 | 1.01 | 1.10 |  |  |
| He | Jagau 2016, FCI | `default` |  |  |  |  | 1.20 | 1.10 | 1.02 | 1.08 | 1.16 | 1.18 |
| Ne | Jagau 2016, CCSD, all electrons | `default` |  |  |  | 0.37 | 0.71 | 0.66 | 0.82 | 0.94 | 0.90 |  |
| Ne | Bahl et al. 2016, one electron | `one_electron` | 1.34 | 1.30 | 1.26 | 1.23 | 1.17 | 1.13 | 1.09 |  |  |  |
| Ar | Jagau 2016, CCSD, all electrons | `default` | 0.17 | 0.59 | 0.78 | 0.99 | 1.15 | 1.16 | 1.11 |  |  |  |
| Ar | Zhang et al. 2014, one electron, Tong-Lin potential | `one_electron` | 0.94 |  |  |  | 0.95 |  | 1.01 |  |  |  |
| Ar | Zhang et al. 2014, one electron, GSZ potential | `one_electron` | 1.12 |  |  |  | 1.16 | 1.20 | 1.17 |  |  |  |
| Ar | Zhang et al. 2014, one electron, Muller potential | `one_electron` | 1.19 |  |  |  | 1.18 |  | 1.16 |  |  |  |
| Ar | tRecX, one electron, Tong-Lin potential | `one_electron` | 0.94 |  |  |  | 0.95 | 1.01 | 1.01 |  |  |  |
| Ar | Bahl et al. 2016, one electron | `one_electron` | 0.97 | 0.93 | 0.91 | 0.91 | 0.95 | 0.99 | 1.01 |  |  |  |
| Kr | Bahl et al. 2016, one electron | `one_electron` | 1.12 | 1.12 | 1.15 | 1.18 | 1.22 | 1.18 |  |  |  |  |
| Kr | tRecX, one electron, Zhang et al.'s GSZ potential | `one_electron` |  | 1.08 |  | 1.15 | 1.17 |  |  |  |  |  |
| Kr | Santra et al. 2006, spin-orbit, all 4p electrons | `default` |  |  |  |  | 1.04 |  |  |  |  |  |
| Xe | Bahl et al. 2016, one electron | `one_electron` | 0.88 | 0.93 | 1.00 | 1.04 | 1.03 |  |  |  |  |  |
| Xe | tRecX, one electron, Zhang et al.'s Tong-Lin-form potential | `one_electron` |  | 0.94 |  | 1.03 |  |  |  |  |  |  |
<!-- ionisation-rate-static:end -->

A column holds the reference's field nearest to its heading, within 7 %. The model
is compared up to twice the barrier-suppression field, the limit of Tong & Lin's
correction. Jagau's rates of neon below 0.14 au and of argon below 0.07 au are those
discussed above.

```{figure} /auto_gallery/ionisation/images/sphx_glr_plot_ionisation_static_001.*
:target: ../auto_gallery/ionisation/plot_ionisation_static.html
:alt: The static rates of helium, neon and argon against their references

The static rates of helium and neon (all electrons) and argon (one electron and all
electrons), and the ratio of each model to the references.
```
