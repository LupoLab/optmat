# Raman lines

Read {doc}`/user-guide/nonlinear-conventions` for the division between electronic
and delayed response. This chapter defines the line parameters; {doc}`raman_response`
assembles them into a causal kernel. Examples start in {doc}`quickstart`, and
{doc}`/user-guide/raman-validation` collects checks and unresolved questions.

The delayed, nuclear part of a molecular gas's third-order response comes from its
Raman lines:
- the rotational S branch, J → J + 2, of a linear molecule;
- the vibrational Q branch of each Raman-active mode.

optmat turns the curated spectroscopy into those lines, each line's frequency, coupling
amplitude and width, at the temperature, density, composition and pump wavelength
asked for (TC-439, TC-442, TC-492, TC-541, TC-542). The response function built from
them, with its tensor form and its density dependence, is on {doc}`raman_response`
(TC-685); its sampling belongs to the propagation code.

```python
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import raman_lines, raman_parameters

lines = raman_lines(
    "N2", 295.0, AMAGAT, {"self": 0.78, "O2": 0.21, "Ar": 0.01}, wavelength=800e-9
)
lines.branch, lines.j  # "S" or "Q", and each line's lower level
lines.frequency  # Omega, rad/s
lines.amplitude  # A, m^5 V^-2 s^-1 (the convention below)
lines.fwhm  # Hz; the dephasing rate is Gamma = pi FWHM
raman_parameters("N2").reasons  # why each constant, strength and width was chosen
```

## The convention

For a field linearly polarised along x, the Raman polarisation of a gas of number
density N is

$$
P_x(t) = \varepsilon_0 N E_x(t) \int_0^\infty h(\tau) E_x^2(t - \tau)\, d\tau,
\qquad h(\tau) = \sum_k A_k \sin(\Omega_k \tau) e^{-\Gamma_k \tau}.
$$

The amplitudes are:

- **Rotational S branch:**
  $A_J = \frac{2}{15} \frac{\Delta\alpha^2}{\hbar\varepsilon_0} \frac{(J+1)(J+2)}{2J+3} (n_J - n_{J+2})$,
  with $n_J$ the population of one sublevel.
  - The nuclear-spin weights follow from exchange symmetry: H2 1:3 (even:odd), D2 and
    N2 6:3, O2 odd N only, CO2 even J only.
  - The anisotropy may depend on J (H2).
- **Vibrational Q branch:**
  $A_J = |\langle 0|\bar\alpha|1\rangle|^2 P_J / (\hbar\varepsilon_0)$, one line per
  rotational level, where $P_J$ is the fraction of molecules in level J.
  - The line positions are $\nu_0 + (B_1 - B_0)J(J+1) - (D_1 - D_0)[J(J+1)]^2$.
  - A band whose rotational structure is not resolved is one line: the spherical tops
    (CH4, SF6), and the bands of CO2 and N2O.
    In `raman.toml`, lower rotational constants imply a resolved band unless
    `resolved = false` is stated. A resolved band needs both ladders, with the upper
    one specified by `upper`, `fit_lines` or `b_difference`; an incomplete request
    raises `DataError` instead of silently producing one unresolved line (TC-735).
  - The transition polarizability comes from one of three sources:
    - a curated transition element, per J where it is tabulated. H2 and D2 use Hunt et
      al.'s ab initio elements, whose mean is (⟨α∥⟩ + 2⟨α⊥⟩)/3.
    - $\partial\bar\alpha/\partial Q$ along the mass-weighted normal coordinate,
      since $|\langle 0|Q|1\rangle|^2 = \hbar/2\omega$. No reduced mass or bond
      convention enters, so polyatomic molecules need no per-bond convention.
    - Schrötter & Klöckner's normalised cross section Σ, the band's cross section
      relative to N2's Q branch with the ν⁴ and Boltzmann factors removed. With N2's
      absolute value $C = 5.05\times10^{-48}$ cm⁶/sr and the band's depolarisation
      ratio ρ, $\bar a^2 = \Sigma C/(16\pi^4)/(1 + 7\rho/(3 - 4\rho))$.
  - **Dispersion:** where Bischel & Black's single-resonance law is curated (N2, O2, H2,
    D2, CH4), the element is scaled from the wavelength at which it holds to the pump,
    by $(\nu_i^2 - \nu_{\rm ref}^2)/(\nu_i^2 - \nu_p^2)$. Pass `wavelength` to
    `raman_lines`. For H2 this law follows Rychlewski's ab initio dispersion of the Q(0)
    element to 0.8 % from 238 to 911 nm (TC-542). The anisotropies are not dispersed.
- **Widths:** $\Delta\nu = c + A(T)/n + n\sum_p x_p\gamma_p(J, T)$, in
  {mod}`optmat.gas.linewidths`:
  - $A$ is the Dicke term of H2 and D2;
  - $\gamma_p$ is the broadening by collision partner p, per line where the data are
    per line, with its temperature law.

  The widths are per line: the overlap and mixing of lines at high density, the
  collapse of a Q branch, is not modelled.

The steady-state Stokes gain of a line follows as $g = \omega_s N A / (2 c^2
\varepsilon_0 \Gamma)$ ({func}`~optmat.gas.raman.stokes_gain`).

## The parameters

Each gas's choices, with the reason for each, are in `data/gas/raman.toml`:

<!-- raman-parameters:start -->
| Gas | B (cm-1) | D (cm-1) | Delta alpha (au) | S FWHM (GHz) | nu_0 (cm-1) | alpha_01 (au) | Q FWHM (GHz) |
|---|---|---|---|---|---|---|---|
| N2 | 1.98957 | 5.760e-06 | 4.521 | 2.46-3.85 | 2329.917 | 0.3789 | 1.36-3.93 |
| O2 | 1.43767 | 4.839e-06 | 6.883 | 1.05-2.95 | 1556.391 | 0.3652 | 1.89-3.21 |
| H2 | 59.33920 | 4.599e-02 | 2.104-2.318 (per J) | 0.08-0.11 | 4160.200 | 0.7392-0.7546 (per J) | 0.378 |
| D2 | 29.91050 | 1.134e-02 | 2.008 | 0.05-0.16 | 2993.500 | 0.6163-0.6229 (per J) | 0.183-0.274 |
| CH4 |  |  |  |  | 2916.700 | 1.1235 | 8.58 |
| SF6 |  |  |  |  | 775.000 | 0.7293 | 48 |
| CO2 | 0.39026 | 1.270e-07 | 14.253 | 5.40-6.60 | 1388.000 | 0.3995 | 6.6 |
| N2O |  |  |  |  | 1285.000 | 0.5404 | 5.57 |
<!-- raman-parameters:end -->

(Widths at 295 K and 1 atm of the pure gas, over the lines kept. The transition
polarizabilities are as curated: static for H2 and D2, near 436 nm for N2 and O2, and
visible for the others.)

- **N2:**
  - Bendtsen's 1974 constants.
  - Wahlstrand et al.'s anisotropy, measured from the rotational response at 800 nm.
  - Murphy et al.'s ∂α/∂Q. Schrötter & Klöckner's absolute cross section gives
    0.375-0.384 au, which brackets it (TC-518).
  - Herring et al.'s per-J rotational widths and temperature laws.
  - **Q widths** (TC-541): Rosasco et al.'s measured Q(J) widths at 298 K, with a
    power law per J fitted to Rahn & Palmer's tables at 500-1500 K. These replace
    Koszykowski et al.'s exponential-gap law, and they reproduce every table to an rms
    of 0.6-5.5 %.
- **O2:**
  - Hajigeorgiou's levels, which Yu et al.'s Dunham fit reproduces (TC-523).
  - Wahlstrand et al.'s anisotropy.
  - Hölzer et al.'s per-N rotational widths, self and in air, from 295 to 1900 K.
  - Millot et al.'s per-N Q widths, self and by N2.
- **H2:**
  - Stoicheff's (1957) B0, D0 and H0, from his paper (TC-586). Weber (1979) reprints
    them, with a misprinted B0 of 56.3392 for 59.3392 (TC-515). These reproduce Teal &
    MacWood's lines to within 0.9 cm⁻¹, their scatter. Komasa et al.'s (2011) ab initio
    lines agree with Stoicheff's within his stated errors.
  - Bishop & Pipin's ab initio anisotropy per J.
  - Hunt et al.'s ab initio static transition elements per J (0.739-0.755 au), with
    Bischel & Black's dispersion.
  - Bischel & Dyer's Q-branch law.
  - Herring et al.'s rotational widths.
- **D2:**
  - Stoicheff's constants.
  - The ab initio optical anisotropy.
  - Hunt et al.'s elements, whose Q(2) agrees with Murphy et al.'s ∂α/∂Q to 0.02 %.
  - **Rotational widths:** Keijser et al.'s (1974) S0(J) widths (their half widths,
    doubled), which Knaap et al.'s figures confirm to 5 % (TC-586). Smyth et al.'s
    reprint of them differs by the rounding of its conversion (at most 0.4 %). Edwards
    et al.'s are 1.3-2.6 times lower, with no single factor
    (TC-530, TC-541).
  - **Q widths:** Smyth et al.'s per-J widths, self and by He, Ar, H2 and CH4, with
    Russell & Roh's Dicke term. Above a few amagat, where the density-dependent term dominates, Q(1) is 1.5 times
    narrower than Q(2) and Q(0) is 1.4 times broader; at 1 amagat both ratios are 1.24.
- **CH4:**
  - Bermejo et al.'s ν1.
  - The mean of Schrötter & Klöckner's thirteen Σ (8.55, ρ = 0), 5 % above Murphy et
    al.'s ∂α/∂Q, with Bischel & Black's dispersion.
  - Ottusch & Rockwell's width law.
- **SF6:**
  - Murphy et al.'s ν1.
  - Schrötter & Klöckner's Σ (3.3-3.96), 3 % above Murphy et al.'s. The four values
    span 20 %, and Everall et al.'s 248 nm cross section needs 23 % more than the
    largest (TC-517).
  - Everall et al.'s width at 5 atm, taken as independent of density.
- **CO2:**
  - Wang's B0 and D0, as Weber reprints them.
  - Bogaard et al.'s optical anisotropy at 632.8 nm, as Gentle et al. reprint it.
  - Pinter's per-line S widths, as Jammu et al. reprint them.
  - The two members of the ν1/2ν2 Fermi dyad (1388 and 1285 cm⁻¹) as two unresolved
    bands. Their positions and ρ are Danichkin's recommended values, and their elements
    come from Schrötter & Klöckner's Σ.
  - **Stopgap:** no Q-branch width is curated, so the S-line coefficient stands in.
- **N2O** (vibrational only):
  - ν1 (1285 cm⁻¹) and ν3 (2224 cm⁻¹) at Danichkin's positions, with elements from
    Schrötter & Klöckner's Σ.
  - **Stopgaps:**
    - ρ = 0 is taken where no depolarisation ratio is curated, which puts an upper
      bound on the element.
    - The width is Chen et al.'s rotational dephasing time (23.8 ps at 2.4 atm), taken
      as collisional.
  - **Not yet:** no rotational lines. No source of N2O's rotational constants is
    curated: the Weber volume has none, and Chen et al.'s B is a simulation input.

<span id="how-well-the-lines-give-the-measured-gains"></span>
<span id="the-rotational-nonlinear-index"></span>
<span id="open-questions"></span>
<span id="against-solprop"></span>

See {doc}`/user-guide/raman-validation` for comparisons with measured quantities.
