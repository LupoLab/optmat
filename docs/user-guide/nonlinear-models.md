# Third-order models and fits

Read {doc}`nonlinear-conventions` first. This chapter explains the polynomial model,
its time-domain form and fit choices. Use {doc}`nonlinear-support` and
{doc}`/nonlinear_accuracy` to assess a calculation. For a runnable starting point,
return to {doc}`/nonlinear`; commands and configuration maintenance are in
{doc}`/developers/fitting`.

## The dispersion model

Below resonance gamma is analytic in the squared frequencies. Bishop showed, and
Shelton's measurements across processes confirm, that to fourth order it depends on the
four frequencies only through $W_2$ and $W_4$:

$$
\gamma = \gamma_0 \left[1 + A\,W_2 + B\,W_2^2 + B'\,W_4\right].
$$

- **$\gamma_0$** is the static gamma, in au.
- **$A$** (in $10^{-10}$ cm$^2$) is the same for every process. To second order in
  frequency one curve against $W_2$ serves all processes.
- **$B$ and $B'$** (in $10^{-20}$ cm$^4$) are the fourth-order terms.
  - $B'$ separates the processes: Bishop & De Kee (1996) give it for helium, and
    calculations of several processes fix it for the other gases.
  - Shelton's single process-dependent $B$ is this $B + B' W_4/W_2^2$. For the ESHG,
    $W_4/W_2^2 = 1/2$; for the Kerr effect, 1/4; for the third harmonic, 7/12.
  - Where no calculation gives the process dependence, $B'$ is held at its prior, 0.

$W_2$ is the natural variable: THG at 1055 nm has the same $W_2$ as the Kerr effect at
609 nm, and so nearly the same gamma. The stored model is valid over its declared
fundamental-wavelength range and, for any process, while $W_2$ stays within the fitted
data. {meth}`~optmat.gas.nonlinear.HyperpolarizabilityModel.gamma_au` warns
({class}`~optmat.exceptions.ModelRangeWarning`) outside the wavelength range, and
wherever the process's $W_2$ exceeds the largest fitted $W_2$ by more than 10 %
(`optmat.gas.nonlinear.W2_MARGIN`): the third harmonic of a fundamental shorter than the
table's THG limit below, for instance, is an extrapolation of the form, which fails
toward resonance.

The fitted parameters (gamma in au; A in 1e-10 cm²; B, B' in 1e-20 cm⁴), with the
largest $W_2$ among the fitted values and the shortest Kerr and THG fundamentals it
corresponds to:

<!-- nonlinear-fits:start -->
| Gas | gamma_0 (au) | A (1e-10 cm²) | B (1e-20 cm⁴) | B' (1e-20 cm⁴) | Values fitted | chi² / dof | Largest W2 fitted (1e9 cm⁻²) | Shortest Kerr / THG fundamental (nm) |
|---|---:|---:|---:|---:|---:|---|---:|---|
| He | 43.15 | 0.448 | 0.210 | -0.043 | 52 | 1.24 / 48 | 5.78 | 263 / 456 |
| Ne | 107.1 | 0.550 | 0.298 | -0.003 | 26 | 7.11 / 21 | 2.86 | 374 / 648 |
| Ar | 1168 | 1.115 | 1.780 | -0.233 | 41 | 24.4 / 36 | 2.86 | 374 / 648 |
| Kr | 2602 | 1.364 | 3.552 | 0.000 | 6 | 20.3 / 2 | 2.52 | 399 / 690 |
| Xe | 6836 | 1.703 | 7.807 | -0.838 | 10 | 12 / 5 | 2.52 | 399 / 690 |
| H2 | 684.9 | 1.288 | 1.899 | -0.151 | 56 | 29.4 / 51 | 2.86 | 374 / 648 |
| D2 | 667.2 | 1.267 | 1.848 | -0.179 | 46 | 13.5 / 41 | 2.64 | 389 / 674 |
| N2 | 920.3 | 1.050 | 1.563 | -0.297 | 37 | 38.7 / 32 | 2.86 | 374 / 648 |
| O2 | 962.8 | 1.112 | 4.602 | 0.000 | 20 | 5.7 / 16 | 2.64 | 389 / 674 |
| CO2 | 1151 | 1.219 | 2.545 | -0.203 | 11 | 15.1 / 6 | 2.52 | 399 / 690 |
| CH4 | 2590 | 1.532 | 4.333 | 0.000 | 13 | 2.04 / 9 | 2.64 | 389 / 674 |
| SF6 | 1429 | 0.494 | 2.154 | 0.000 | 8 | 3.37 / 4 | 2.86 | 374 / 648 |
| CO | 1432 | 1.402 | 1.489 | 0.000 | 6 | 0.638 / 2 | 2.49 | 401 / 695 |
| N2O | 1804 | 1.250 | 2.400 | 0.000 | 1 | 0 / -2 | 0.907 | 664 / 1150 |
| O3 | 2712 | 0.728 | 3.544 | 0.000 | 4 | 0 / 1 | 0.829 | 695 / 1203 |
<!-- nonlinear-fits:end -->

### The time-domain form

Because gamma depends on each frequency only through its square and the square of the
sum, the model has an exact time-domain form. The operator

$$
\hat D = -\frac{1}{(2\pi c)^2}\frac{\partial^2}{\partial t^2}
$$

multiplies a field component of wavenumber $\nu$ by $\nu^2$. Applied to one input field
it gives $\nu_i^2$; applied to the product it gives $\nu_\sigma^2$. The polarization of
the model, for any real field $E(t)$, is then

$$
\begin{aligned}
P = \frac{N\gamma_0}{6}\Big\{ & E^3
 + A\left[\hat D(E^3) + 3E^2\,\hat DE\right] \\
 & + B\left[\hat D^2(E^3) + 3E^2\,\hat D^2E + 6\,\hat D\!\left(E^2\,\hat DE\right)
 + 6E\,(\hat DE)^2\right]
 + B'\left[\hat D^2(E^3) + 3E^2\,\hat D^2E\right]\Big\},
\end{aligned}
$$

with $A$, $B$, $B'$ in consistent units. Each term reproduces, for every mixing channel
at once (Kerr, cross-phase, the third harmonic and every four-wave-mixing product), the
gamma of that channel's frequencies. Every term is a product of time derivatives of the
field and of field products, so on a grid it costs a few FFTs. The polynomial grows as
$\nu^4$: it must only see the band where the model is valid, so the field should be
limited to that band or the terms cut off above it (*Where the fits are supported*,
in {doc}`/user-guide/nonlinear-support`). A propagation code can use this form
subject to the same support and resonance limitations.

## The fits

<span id="configuration"></span>

Configuration fields are documented in {doc}`/developers/fitting`.

### The fit

{func}`~optmat.gas.nonlinear.fit_hyperpolarizability` converts every fitted dataset to
electronic $\gamma_T$ and its frequency sums, then minimises

$$
\chi^2 = \sum_{d}\sum_{i\in d}\left(\frac{s_{g(d)}\,\gamma(W_{2,i}, W_{4,i}) -
\gamma_i}{\sigma_i}\right)^2 + \sum_\text{priors}\left(\frac{p - p_0}{\sigma_p}\right)^2
$$

by nonlinear least squares:
- **Scales:** each dataset with `free_scale` or `systematic` belongs to a scale group
  $g$, with scale $s_g$ (1 for the rest). A systematic prior enters as $(s_g - 1)/
  \sigma$.
- **Priors:** on $B'$, and where given on $A$ and $B$.
- **Covariance:** the inverse of $J^TJ$ at the solution, over the four physical
  parameters.
- **Statistics:** $\chi^2$ and the degrees of freedom (values minus physical parameters
  minus scales; the priors are not counted).
- **Uncertainty:** {meth}`~optmat.gas.nonlinear.HyperpolarizabilityModel.gamma_sigma_au`
  propagates the covariance to any process and wavelength. It scales it by the reduced
  $\chi^2$ where that exceeds one (Birge).

The data that are not fitted are validators: {func}`~optmat.gas.nonlinear.validation`
reports, per validator, the mean ratio of data to model and its $\chi^2$, or why it
could not be compared (a CARS value without its Stokes wavelength).

<span id="running-and-storing"></span>

Commands, stored records and regeneration order are in {doc}`/developers/fitting`.

### The choices, gas by gas

- **The backbone:** Shelton's 1990 ESHG, on Bishop & Pipin's helium scale. It gives his
  own 1064 and 1319 nm values and his rescaling of the visible literature to that scale
  (`shelton_1990/table3_<gas>_present` and `_literature`, and Shelton 1992 for neon).
- **Calculated process dependence,** one free scale per calculation (`scale_group`), so
  that it informs $B'$ and the shape but not the level:
  - Hättig & Jørgensen's CCSD Kerr, ESHG and THG for Ne, Ar and N2;
  - Bishop et al. 1991 (electronic) for H2 and D2;
  - Tarazkar et al. (the q-aug-cc-pVQZ-pp rows) for Xe;
  - Sekino & Bartlett's CCSD(T) for CO2 and CO.

  A calculated THG enters only for fundamentals of 650 nm and longer, where the
  fourth-order form holds. The calculations' nearest fundamental above that is 694.3 nm
  for all five gases, so the harmonics that enter are at 231.4 nm and longer (not 217 nm,
  as the manual said before; TC-864).
- **He:** on theory alone: Bishop & Pipin's four processes to 0.1 hartree (455.6 nm),
  with 0.1 % allowed per value. The measurements are validators.
- **Kr, O2, CH4, SF6:** Shelton only, with $B'$ held at its prior. CH4's 1319 nm point
  is excluded, as Shelton excluded it. Fitting Hättig & Jørgensen's CCSD shape for CH4
  was tried and rejected (scale 0.87, $\chi^2$ 350 for 18; TC-314).
- **CO:** Ward & Miller's ESHG at 694 nm sets the level; Sekino & Bartlett's four
  processes the dispersion.
- **N2O:** only Wahlstrand et al.'s n2 at 800 nm, so the dispersion is CO2's, as priors
  of 50 % on $A$ and 100 % on $B$.
- **O3:** theory only: Naves et al.'s CCSD(T) static and Kerr values, with 12 % allowed
  for the method dependence, from 700 nm. Its third harmonic, for fundamentals up to
  about 1050 nm, falls in the Hartley and Huggins bands and is not given.

The gallery page *The second hyperpolarizability: every process on one curve* shows the
data of every process against $W_2$ with the fitted curve.

## Values

At 800 nm and 1800 nm, with the statistical one-sigma uncertainty of the 800 nm values.
The working uncertainties, which add the leave-one-out spread and the systematic floors,
are on the accuracy page.

<!-- nonlinear-values:start -->
| Gas | gamma_0 (au) | Kerr, 800 nm (au) | THG, 800 nm (au) | Kerr, 1800 nm (au) | n2, 800 nm, 1 atm, 20 °C (1e-20 cm²/W) |
|---|---:|---:|---:|---:|---:|
| He | 43.15 | 44.39 ± 0.0069 | 47.06 ± 0.017 | 43.39 | 0.368 |
| Ne | 107.1 | 110.9 ± 0.26 | 119.3 ± 0.27 | 107.8 | 0.92 |
| Ar | 1168 | 1258 ± 1.9 | 1480 ± 2.3 | 1185 | 10.4 |
| Kr | 2602 | 2860 ± 17 | 3592 ± 31 | 2647 | 23.8 |
| Xe | 6836 | 7767 ± 30 | 1.078e+04 ± 77 | 6988 | 64.8 |
| H2 | 684.9 | 745 ± 0.75 | 893.9 ± 1.2 | 696 | 6.18 |
| D2 | 667.2 | 724.7 ± 0.62 | 866.6 ± 1 | 677.8 | 6.01 |
| N2 | 920.3 | 986.1 ± 1.7 | 1147 ± 1.6 | 932.5 | 8.19 |
| O2 | 962.8 | 1047 ± 2.5 | 1319 ± 5.1 | 976.7 | 8.7 |
| CO2 | 1151 | 1249 ± 6.3 | 1512 ± 5.4 | 1168 | 10.4 |
| CH4 | 2590 | 2882 ± 8.7 | 3729 ± 11 | 2641 | 24 |
| SF6 | 1429 | 1485 ± 4.6 | 1670 ± 4.5 | 1438 | 12.5 |
| CO | 1432 | 1566 ± 44 | 1884 ± 53 | 1458 | 13 |
| N2O | 1804 | 1961 ± 2.5e+02 | 2378 ± 3.4e+02 | 1832 | 16.4 |
| O3 | 2712 | 2873 ± 3e+02 | n/a (resonant) | 2738 | 23.9 |
<!-- nonlinear-values:end -->
