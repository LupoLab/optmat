# Ionisation rates of molecules

The rate at which an intense, linearly polarised field ionises a molecule of an
**unaligned gas**: {mod}`optmat.gas.molecular_ionisation_rates`. It is the rate of the
noble gases ({doc}`ionisation_rates`) carried to molecules, with a correction at
finite field fitted to all-electron static rates. Beside it are the forms of the
literature, which propagation codes use, for comparison and to reproduce results
computed with them ({mod}`optmat.gas.molecular_ionisation_validation`).

```python
import numpy as np

from optmat.gas.molecular_ionisation_rates import molecular_rate_model

field = 2.74e10 * np.cos(np.linspace(0.0, 2 * np.pi, 101))  # one cycle, V/m

nitrogen = molecular_rate_model("N2", 800e-9)  # the default
nitrogen.cycle_averaged_rate(2.74e10)  # 1/s, at the amplitude of 1e14 W/cm^2
nitrogen.instantaneous_rate(field)  # 1/s, along a carrier-resolved field in V/m
nitrogen.electrons, nitrogen.tong_lin_alpha, nitrogen.structure_factors
nitrogen.ionization_energy  # J: the vertical energy, which the rate uses
nitrogen.adiabatic_ionization_energy  # J: the energy of the ion's ground state
nitrogen.citation, nitrogen.sources, nitrogen.reasons, nitrogen.uncertainty_factor
molecular_rate_model("N2", 800e-9, "moadk_zhao_2010")  # a form of the literature
```

The molecules are H2, D2, N2, CO, O2, CO2, N2O, CH4 and SF6. The parameters are those
of {doc}`ionisation`. The field is in V/m, the rate in s⁻¹.

## The rate

Atomic units: $\kappa = \sqrt{2 I_p}$, $F_0 = \kappa^3$, $n^* = 1/\kappa$, the
Keldysh parameter $\gamma = \omega\kappa/F$.

$$
W(F, \omega) = N_\text{eff}\, V(F) \sum_{m'} K_{m'}\, w^\text{sr}_{m'}(F, \omega)\,
  Q(F)\, S(F),
\qquad
K_{m'} = \frac{\langle B^2(m') \rangle\, \kappa^{1 - 2n^*}}{2^{|m'|}\,|m'|!\, I_p} .
$$

| Factor | What it is | From | Register |
|---|---|---|---|
| $\langle B^2(m') \rangle$ | the structure factor of the ionising orbital for the projection $m'$ on the field, averaged over the orientations of the molecule | Tong, Zhao & Lin 2002, eqs. (8)-(9); coefficients of Zhao et al. 2010 (2017 for N2O, 2011 for CH4) | TC-673, TC-677 |
| $w^\text{sr}_{m'}$ | the rate of Perelomov, Popov & Terent'ev in a short-range potential for the projection $m'$, with $K_{m'}$ in the place of an atom's $\lvert C_{\kappa l}\rvert^2 f(l, m)$ | Perelomov et al. 1966, eq. (54) | TC-677 |
| $Q$ | the static Coulomb factor $(2F_0/F)^{2n^*}$ | Perelomov & Popov 1967, eq. (42) | TC-653, TC-670 |
| $S$ | $\exp[-\alpha F/(I_p \kappa^3)]$, Tong & Lin's form of the correction at finite field | $\alpha$ fitted here | TC-678 |
| $N_\text{eff}$ | the number of electrons in the ionising shell | the weak-field limit (Tolstikhin et al. 2014) | TC-678 |
| $V$ | the average of the rate over the vibrational ground state, $e^{c h(\gamma)^2/F^2}$ | $c$ fitted here to Saenz 2002 (H2) | TC-679 |

- **For one partial wave** $K_{m'}$ is an atom's $|C_{\kappa l}|^2 f(l, m')$ averaged
  over $m$, and the rate is the noble gases' rate: the tests hold the two to 1e-9 for
  helium's constants.
- **Orientation.** The rate is the average over orientations of the rate. For N2 the
  ionised fraction of a 10 or 40 fs pulse it gives is within 2 % of the average over
  orientations of the ionised fractions, at any depletion, with the all-electron
  rates' dependence on the angle; the alignment the pulse itself induces changes the
  rate by less than 7 % for 5 to 50 fs (TC-676). What no rate $W(|E|)$ holds is the
  alignment a pulse leaves behind: a second pulse within a rotational revival (8.4 ps
  for N2, 280 fs for H2) meets a gas whose rate differs by up to the anisotropy of the
  table below.
- **The energy** is the vertical one, that of the neutral's geometry, which a
  tunnelling electron at fixed nuclei sees. The adiabatic energy is kept as
  `adiabatic_ionization_energy`: what an ionisation costs the field is for the
  propagation code to decide (TC-679).
- **The nuclei** are fixed at the neutral's geometry, except for H2 and D2, whose rate
  is multiplied by $V$. For N2, Kornev & Zon's (2015) sum over the ion's vibrational
  levels is the fixed-nuclei rate at the vertical energy to 2 % (TC-682).
- **The two rates.** The cycle-averaged and the instantaneous rate are related as for
  the noble gases, through the rate's own average over a cycle (TC-662, TC-828).
- **Above the correction's range** (each molecule's `correction_range`, 1.61 to 1.96
  times its barrier-suppression field: the largest field its alpha was fitted to or
  checked against) the rate is held at its value there, with a warning (TC-828).

<!-- molecular-rate-constants:start -->
| Molecule | N_eff | alpha | I_p, rate (eV) | I_p, adiabatic (eV) | Coefficients | c (au) | Factor |
|---|---|---|---|---|---|---|---|
| H2 | 2 | 6 | 16.060 | 15.426 | zhao_2010 | 0.00165 | 1.5 |
| D2 | 2 | 6 | 16.060 | 15.466 | zhao_2010 | 0.00165 | 2 |
| N2 | 2 | 14 | 15.580 | 15.581 | zhao_2010 |  | 1.5 |
| CO | 2 | 14 | 14.010 | 14.010 | zhao_2010 |  | 1.5 |
| O2 | 2 | 6.1 | 12.307 | 12.070 | zhao_2010 |  | 2 |
| CO2 | 4 | 14 | 13.790 | 13.778 | zhao_2010 |  | 3 |
| N2O | 4 | 14 | 12.890 | 12.886 | zhao_2017 |  | 3 |
| CH4 | 6 | 18 | 14.400 | 12.610 | zhao_2011 |  | 2 |
| SF6 | 6 | 14 | 15.690 | 15.350 | zhao_2011 |  | 10 |
<!-- molecular-rate-constants:end -->

"Factor" is the factor within which the rate is expected to be right, with its basis
in `reasons`. The constants:

- **$N_\text{eff}$** is the number of electrons in the ionising shell, each orbital of
  a degenerate shell counted: 2 for the $\sigma_g$ shells of H2 and N2 and the
  half-filled $\pi_g$ of O2, 4 for the $\pi$ shells of CO2 and N2O, 6 for the $t_2$
  shell of CH4. It is the weak-field limit of the many-electron theory, and the count
  Kjeldsen et al. (2005) apply. Fitted freely, the all-electron rates ask for 2.1 (H2)
  and 1.8 (N2) with Zhao et al.'s coefficients of 2010, and for 2.3 and 1.1 with those
  of 2017: the coefficients and the count are only determined together, and the
  choice of 2010 is the one that makes the count the theory's (TC-678).
- **$\alpha$** is fitted with that count to the average over orientations of Majety &
  Scrinzi's (2015) all-electron static rates at fixed nuclei: 6.0 for H2 (Tong & Lin's
  own value for H2⁺ is 6) and 14.0 for N2; and to the average of Jagau's (2018)
  coupled-cluster rates: 6.1 for O2 and 18.0 for CH4 (TC-698). CO, CO2, N2O and SF6
  take N2's. The constant is not one number: O2's $\pi_g$ orbital asks for H2's 6,
  not N2's 14, and whether CO2's and N2O's $\pi$ orbitals follow O2 is not known. For CO,
  isoelectronic with N2, all-electron rates exist and test the transfer (below); for
  the others none is published.
- **$V$**: the bond length enters the rate through the ionisation energy, and the
  rate depends on the energy chiefly through its exponent $(2F_0/3F)\,g(\gamma)$, whose
  derivative is $-(2\kappa/F)\,h(\gamma)$, $h = g + \gamma g'/3$. For a Gaussian
  spread of the energy over the ground state the average of the rate over it is
  exactly $e^{c h^2/F^2}$: $h = 1$ in a static field, and $h \approx \ln(2\gamma)/\gamma$
  in the multiphoton limit, where the rate depends on the energy only through the
  number of photons. The average is dominated by the bond length $x = \sqrt{2c}\,h/F$
  standard deviations out in the ground state; beyond two, where neither the
  Gaussian nor the linear dependence of the energy holds, $\ln V$ is continued as
  $2x - 2$, the contribution of that bond length, so that the rate still falls to
  zero with the field (TC-679).
- **$c$** is fitted to Saenz's (2002) rates of H2 averaged over its vibrational ground
  state, over his rates at 1.4 bohr, from 0.06 to 0.12 au: 0.00165. D2 takes H2's, an
  assumption: in the harmonic estimate its $c$ is smaller by $\sqrt 2$.
- **CH4's energy** is the middle maximum of a band split by the Jahn-Teller effect
  (13.6, 14.4 and 15.0 eV): the vertical energy of the degenerate state at the
  neutral's geometry, which Jagau's calculation puts at 14.43 eV. With the first
  maximum no $\alpha$ reproduces the dependence of his rates on the field (TC-698).

**SF6's rate is a prediction that nothing tests.** Its coefficients (Zhao et al.
2011) are read on the page, as the paper library's extraction of that part of their
table is unusable (TC-684). They are dominated by an l = 4 wave: the HOMO is a
combination of fluorine lone pairs 3 bohr from the centre, and a one-centre expansion
of such an orbital has large high partial waves. The orientation-averaged structure
factor at m′ = 0 is 75 per orbital, against 4.2 for N2 and 5.7 for CH4. With 6
electrons against N2's 2, the default is 80 times N2's rate at 2e14 W/cm² for nearly
the same ionisation energy, and 100 to 2300
times the hydrogen-like `adk` form. No calculation or measurement of SF6's
rate is in the library to decide it; its factor is 10.

**The range.** The rate is given without a warning up to 1.5 times the
barrier-suppression field ($I_p^2/4$; `maximum_field`), for 350 nm and longer
wavelengths. To 1.5 times, N2's correction holds its all-electron rates within 13 %.
Beyond, the correction outruns the rate: with N2's $\alpha$ of 14 the rate of every
molecule that takes it stops rising at 1.7 times the barrier-suppression field, which
no rate can do. The atoms' smaller constants keep their rates rising to twice it, the
atoms' range (TC-683).

```{figure} auto_gallery/ionisation/images/sphx_glr_plot_molecular_ionisation_default_001.*
:target: auto_gallery/ionisation/plot_molecular_ionisation_default.html
:alt: The default rate of each molecule at 800 nm, N2 at several wavelengths, and H2's vibrational factor

The default at 800 nm for each molecule (circles: Keldysh parameter 1; squares: the
barrier-suppression intensity; band: N2's uncertainty factor); N2 from 400 nm to a
static field; H2's vibrational factor, static and at 800 and 1800 nm, against Saenz's.
```

## The forms of the literature

Atomic units: $\kappa = \sqrt{2 I_p}$, $F_0 = \kappa^3$, $n^* = Z/\kappa$.

<!-- molecular-rate-models:start -->
| Model | Form | Frequency | Molecules |
|---|---|---|---|
| `moadk_zhao_2010` | moadk | static | H2, D2, N2, CO, O2, CO2 |
| `moadk_tong_2002` | moadk | static | H2, D2, N2, CO, O2 |
| `moadk_zhao_2017` | moadk | static | H2, D2, N2, CO, O2, CO2, N2O |
| `ppt_talebpour` | ppt | static Coulomb factor | D2, N2, O2 |
| `ppt_talebpour_ilkov` | ppt | ilkov Coulomb factor | N2, O2 |
| `adk_hydrogenic` | adk | static | H2, D2, N2, CO, O2, CO2, N2O, CH4, SF6 |
| `default` | molecular | static Coulomb factor | H2, D2, N2, CO, O2, CO2, N2O, CH4, SF6 |
<!-- molecular-rate-models:end -->

**`moadk`: the molecular tunnelling rate** of Tong, Zhao & Lin (2002). The ionising
orbital enters through the coefficients of its asymptotic tail, $\Psi^m \to \sum_l
C_l\, r^{Z/\kappa - 1} e^{-\kappa r}\, Y_{lm}$ (their eqs. 1-2). For a molecule at the
Euler angles $\mathbf R$ to a static field (eqs. 8-9)

$$
w(F, \mathbf R) = \sum_{m'} \frac{B^2(m')}{2^{|m'|}\,|m'|!}\,
  \frac{1}{\kappa^{2Z/\kappa - 1}}
  \left(\frac{2\kappa^3}{F}\right)^{2Z/\kappa - |m'| - 1} e^{-2\kappa^3/3F},
\qquad
B(m') = \sum_l C_l\, D^l_{m'm}(\mathbf R)\, Q(l, m') ,
$$

with $Q(l, m) = (-1)^m \sqrt{(2l+1)(l+|m|)!/2(l-|m|)!}$. The models here are its
**average over the orientations** of the molecule: the cross terms between partial
waves vanish by the orthogonality of the rotation matrices, and

$$
\langle B^2(m') \rangle = \sum_l C_l^2\, \frac{Q^2(l, m')}{2l + 1}
  = \sum_l C_l^2\, \frac{(l + |m'|)!}{2\,(l - |m'|)!} ,
$$

whatever the orbital's own $m$. The sum over $m'$ runs over both signs. It is a sum
over the projections of one orbital on the field, not a count of electrons.

**`ppt`: an effective charge.** Talebpour, Yang & Chin (1999) fitted measured ion
yields of N2 and O2, and Talebpour, Larochelle & Chin (1998) modelled those of D2,
with the rate of Perelomov, Popov & Terent'ev for an electron with $l = m = 0$ in the
field of a charge $Z_\text{eff}$ below 1: $n^* = Z_\text{eff}/\kappa$ and the Hartree
coefficient with $l^* = n^* - 1$. The charge is a constant of that rate, not a
property of the molecule; with another Coulomb factor or another ionisation energy it
is not the fitted constant.

- **`ppt_talebpour` is the rate they used**, with the static Coulomb factor
  $(2F_0/F)^{2n^*}$ (their rate is Larochelle et al.'s eq. 10).
- **`ppt_talebpour_ilkov` uses the same charges in Ilkov et
  al.'s form**, which is lower by $(1 + \gamma^2)^{n^*}$, 1.8 at $\gamma = 1$ for N2.
  The charges were not fitted within it.
- O2's 0.53 was fitted with $I_p = 12.55$ eV, and both models take 12.07 eV (TC-542).

**`adk`: the ionisation energy alone.** The tunnelling rate of Ammosov, Delone &
Krainov for a hydrogen-like s electron, $Z = 1$, $l^* = 0$, no angular factor:

$$
w = |C_{n^* 0}|^2\, I_p \left(\frac{2F_0}{F}\right)^{2n^* - 1} e^{-2F_0/3F},
\qquad
|C_{n^* 0}|^2 = \frac{2^{2n^*}}{n^*\,\Gamma(n^* + 1)\,\Gamma(n^*)} .
$$

It is what can be written for a molecule of which nothing but the ionisation energy
is known, and it knows nothing of the orbital: for O2, whose $\pi_g$ orbital has a
node along the field, it is 7 to 19 times the other forms.

**The cycle average.** The `moadk` and `adk` forms are static. Their cycle-averaged
rate is the static rate at the amplitude times $(3F/\pi F_0)^{1/2}$, the average in the
tunnelling limit (Perelomov et al. 1966, eq. 7). The `ppt` form is cycle averaged, and
its instantaneous rate is that divided by the same factor. For the noble gases optmat
replaces this factor by one that holds at any Keldysh parameter (TC-662); the forms
here keep it, as the codes do.

**One electron.** Each form is the rate of one electron, as its source gives it.
To apply an occupancy of two, use `dataclasses.replace(model, electrons=2.0)`.
Keep the occupancy consistent with the coefficients and the quantity being compared.

## The models against one another

The cycle-averaged rate at 800 nm, in s⁻¹: the default's, and the forms' of one
electron:

<!-- molecular-rate-forms:start -->
| Molecule | I (W/cm2) | `moadk_zhao_2010` | `moadk_tong_2002` | `moadk_zhao_2017` | `ppt_talebpour` | `ppt_talebpour_ilkov` | `adk_hydrogenic` | `default` |
|---|---|---|---|---|---|---|---|---|
| H2 | 5e+13 | 2.28e+08 | 4.51e+08 | 2.13e+08 |  |  | 3.31e+08 | 7.37e+09 |
| H2 | 1e+14 | 1.03e+11 | 2.04e+11 | 9.63e+10 |  |  | 1.50e+11 | 3.98e+11 |
| H2 | 2e+14 | 7.48e+12 | 1.48e+13 | 6.99e+12 |  |  | 1.09e+13 | 9.99e+12 |
| D2 | 5e+13 | 2.07e+08 | 4.11e+08 | 1.94e+08 | 1.27e+09 |  | 3.02e+08 | 7.37e+09 |
| D2 | 1e+14 | 9.62e+10 | 1.91e+11 | 8.99e+10 | 9.92e+10 |  | 1.40e+11 | 3.98e+11 |
| D2 | 2e+14 | 7.11e+12 | 1.41e+13 | 6.65e+12 | 3.53e+12 |  | 1.04e+13 | 9.99e+12 |
| N2 | 5e+13 | 4.26e+08 | 2.37e+08 | 6.82e+08 | 4.22e+09 | 1.44e+09 | 2.34e+08 | 1.46e+10 |
| N2 | 1e+14 | 2.13e+11 | 1.19e+11 | 3.42e+11 | 3.06e+11 | 1.52e+11 | 1.17e+11 | 7.38e+11 |
| N2 | 2e+14 | 1.67e+13 | 9.33e+12 | 2.68e+13 | 1.00e+13 | 6.58e+12 | 9.07e+12 | 1.51e+13 |
| CO | 5e+13 | 1.60e+10 | 4.92e+09 | 2.59e+10 |  |  | 7.47e+09 | 2.18e+11 |
| CO | 1e+14 | 3.05e+12 | 9.34e+11 | 4.94e+12 |  |  | 1.41e+12 | 5.63e+12 |
| CO | 2e+14 | 1.20e+14 | 3.65e+13 | 1.95e+14 |  |  | 5.48e+13 | 6.18e+13 |
| O2 | 5e+13 | 3.91e+10 | 5.56e+10 | 5.93e+10 | 8.17e+10 | 4.39e+10 | 4.33e+11 | 2.44e+11 |
| O2 | 1e+14 | 2.50e+12 | 3.55e+12 | 3.79e+12 | 2.12e+12 | 1.43e+12 | 2.63e+13 | 4.04e+12 |
| O2 | 2e+14 | 4.59e+13 | 6.52e+13 | 6.96e+13 | 2.98e+13 | 2.37e+13 | 4.50e+14 | 3.35e+13 |
| CO2 | 5e+13 | 1.38e+10 |  | 1.07e+10 |  |  | 1.23e+10 | 3.28e+11 |
| CO2 | 1e+14 | 2.37e+12 |  | 1.84e+12 |  |  | 2.02e+12 | 7.91e+12 |
| CO2 | 2e+14 | 8.88e+13 |  | 6.87e+13 |  |  | 7.10e+13 | 8.28e+13 |
| N2O | 5e+13 |  |  | 8.02e+10 |  |  | 8.09e+10 | 1.14e+12 |
| N2O | 1e+14 |  |  | 8.23e+12 |  |  | 7.86e+12 | 1.85e+13 |
| N2O | 2e+14 |  |  | 2.13e+14 |  |  | 1.89e+14 | 1.33e+14 |
| CH4 | 5e+13 |  |  |  |  |  | 1.43e+11 | 3.56e+11 |
| CH4 | 1e+14 |  |  |  |  |  | 1.19e+13 | 9.90e+12 |
| CH4 | 2e+14 |  |  |  |  |  | 2.54e+14 | 1.12e+14 |
| SF6 | 5e+13 |  |  |  |  |  | 3.93e+08 | 9.07e+11 |
| SF6 | 1e+14 |  |  |  |  |  | 1.69e+11 | 5.22e+13 |
| SF6 | 2e+14 |  |  |  |  |  | 1.19e+13 | 1.25e+15 |
<!-- molecular-rate-forms:end -->

- **At 5e13 W/cm² the default is above every form**: for N2 by 3.5 (the effective-charge
  form within which Talebpour et al. fitted) to 60 (the hydrogen-like form). The
  Keldysh parameter is 1.6 there, and the default carries the rate's dependence on
  the frequency with the static Coulomb factor, which the static forms lack. At
  2e14 W/cm² the models of N2 come within a factor of 2.3 of one another.
- **For O2 the default is within a factor of 2.2 of the forms that know the orbital**
  from 1e14 to 2e14 W/cm², and above them at 5e13 W/cm², where the frequency counts.
- For **O2** the forms that know the orbital or were fitted to yields agree with one
  another within a factor of 3; the hydrogen-like form is 5 to 19 times above them.
  This is the "suppressed ionisation" of O2 against xenon, an atom of the same
  ionisation energy.
- For **H2 and D2** the tabulations of the structure coefficients differ by a factor
  of 2 in the rate; for **N2** the two of Zhao et al., of 2010 and 2017, by 1.6.

<span id="validation"></span>
<span id="against-all-electron-static-rates"></span>
<span id="against-measured-yields"></span>
<span id="open-questions"></span>

See {doc}`/user-guide/molecular-ionisation-validation`.

## The rates with a dependence on the frequency in the literature

Three forms of a "molecular PPT" rate are printed, and they are not the same (TC-674):

| Source | Coulomb factor on $(2\kappa^3/F)^{2n^* - \lvert m\rvert - 1}$ | Cycle average | Sum over $m'$, orientation |
|---|---|---|---|
| Benis et al. 2004, eq. (2); Zhao, Liu & Zhou 2014, eq. (13) as corrected in 2015 | $(1 + \gamma^2)^{-(2n^* - \lvert m\rvert - 1)/2}$ | none printed | one $m$; the axis along the field |
| Zhao, Liu & Zhou 2014, eq. (13) as printed | $(1 + \gamma^2)^{-(2n^* - \lvert m\rvert - 1)}$: a misprint, by the authors' corrigendum | none printed | the same |
| Sapkota et al. 2026, eq. (5), after Zhao et al. 2016 | $(1 + \gamma^2)^{\lvert m'\rvert/2 + 3/4}\,(1 + 2\gamma/e)^{-2n^*}$ | $(3F/\pi\kappa^3)^{1/2}$ | a sum over $m'$ at one angle, that of the mean $\cos^2\theta$ |

- The first is the form of Ilkov et al. carried to molecules, with the factor $(1 +
  \gamma^2)^{-n^*}$ that Perelomov & Popov's Coulomb correction does not have
  (TC-653, TC-670). The third carries Popruzhenko et al.'s factor.
- None is averaged over the orientations of an unaligned gas: the third evaluates the
  rate at one angle, 54.7° for a gas without alignment, which is the rate at the mean
  angle and not the mean of the rate.
- None has a correction at finite field, a count of electrons, or a statement of
  which ionisation energy is meant.
- None of the papers tabulates a rate or a yield. Their comparisons with measurements
  are of yields normalised to the measurement, which test the dependence on intensity
  and not the magnitude.

None of these is a model here. The default is the rate of {doc}`ionisation_rates`
carried to molecules: the short-range rate of Perelomov et al. with the structure
factors of the molecular orbital averaged over orientations, the static Coulomb
factor, a correction at finite field and the electron count. It is in none of the
papers.

<span id="against-solprop"></span>

See {doc}`/user-guide/molecular-ionisation-validation`.
