# Fitting dispersion models

The package fits its own dispersion models to the curated data, with a physical form,
stated priors and an uncertainty that follows from the data. Read {doc}`models`
and {doc}`conventions` first. This chapter explains the method; {doc}`accuracy`
interprets its uncertainties, and {doc}`/developers/fitting` describes how to run it.

## The model

Below its first resonance the polarizability of an atom is a sum over its dipole
oscillator-strength distribution. The fitted model keeps a few terms of it:

$$
\alpha(\omega) = \sum_{\text{lines}} \frac{f_\ell}{E_\ell^2 - \omega^2}
+ \sum_{j=1}^{M} \frac{F_j}{\mathcal{E}_j^2 - \omega^2}
\qquad\text{(atomic units)},
$$

where the resonance lines sit at their spectroscopic energies $E_\ell$ (NIST ASD, known to
$10^{-7}$) with oscillator strengths $f_\ell$ fitted under priors from the NIST values and
their accuracy grades, and $M$ effective oscillators, of strengths $F_j$ and energies
$\mathcal{E}_j$ above a stated minimum, stand for the rest of the distribution: the higher
Rydberg lines and the continuum. The poles are therefore where the atom's are, the model is
finite and monotonic between the static limit and the first line, and its Cauchy moments
$S(-2k) = \sum f E^{-2k}$ are those of a physical distribution. Priors on dipole sums (the
static polarizability $S(-2)$, the Thomas-Reiche-Kuhn sum $S(0)$) can be added.

## The fit

A fit configuration, `data/gas/fits/<name>.toml`, lists the model's lines and effective
oscillators, the datasets fitted and the options of each, and the datasets used only to
validate. The fit minimises

$$
\chi^2 = \sum_{d}\sum_{i \in d}
\left(\frac{s_d\,\alpha(\omega_i) - \alpha_i}{\sigma_i}\right)^2
+ \sum_\ell \left(\frac{f_\ell - f_\ell^{\text{prior}}}{\sigma_\ell}\right)^2
+ \sum_{d} \left(\frac{s_d - 1}{\sigma^{\text{sys}}_d}\right)^2 ,
$$

in the polarizability, where every value $\alpha_i$ is its dataset's refractivity or
polarizability converted with the density its source meant. Each dataset has a scale $s_d$:
fixed at 1, constrained by the dataset's stated systematic uncertainty, or free. A free scale
takes a dataset's shape and not its level: the choice made for Mansfield & Peck's
helium, whose absolute value is 0.25 % above modern theory (TC-44).

The parameter covariance is $(J^\mathsf{T}J)^{-1}$ at the minimum, and the uncertainty band
of $\alpha(\omega)$ follows by linear propagation, scaled by the reduced $\chi^2$ where that
exceeds one. Several starting points guard against local minima (six by default,
`[model] starts` in the configuration: deuterium needs twenty).

## Other observables

Besides polarizabilities, a fit takes:

- **Cauchy moments** $S(k)$ from ab initio theory (Puchalski et al. for helium, Lesiuk et
  al. for neon and argon), compared with the model's sums directly.
- **Phase-matching density ratios** (Shelton 1986, 1990): two gases phase-match the
  field-induced second harmonic at densities in the ratio
  $\rho_\text{gas}/\rho_\text{ref} = \Delta\alpha_\text{ref}/\Delta\alpha_\text{gas}$,
  $\Delta\alpha = \alpha(2\omega) - \alpha(\omega)$. The partner gas is held at its
  recommended model, so the fits run in dependency order (helium and argon, then neon,
  nitrogen, hydrogen, krypton and xenon, then oxygen and deuterium), and a refit of a
  partner calls for a refit of the gases that use it. A ratio whose partner model is less
  certain than the ratio can bias the fit; nitrogen's ratio to argon is therefore a
  validator (TC-96). `optmat.gas.crosschecks.phase_match_ratios` compares every ratio with the
  recommended models.

Infrared-active molecules carry fixed infrared terms outside the parameter
covariance. Measured refractivities and dielectric constants include the vibrational
response; electronic theory excludes it, and Cauchy moments are electronic.
The CO, CO2, CH4 and N2O fits use the derived band rules and qualified windows
(TC-797). Pure rotation is a separate term class: it enters refractive observations
and is excluded from the fitter's electric observations and moments, as explained
in {doc}`conventions`.

The model exposes the qualified intervals and gap reasons, and
`representation_bound_au(omega, order)` gives its infrared representation-error
bound for orders 0, 2 and 3 inside qualified intervals, with NaN outside. The
{doc}`curated data <data>` page describes the derived tables; {doc}`/developers/fitting` gives
the configuration fields and provenance hashes.

**Rayleigh scattering** cross sections convert to the polarizability's magnitude through
$\sigma = (8\pi/3)\,k^4 (\alpha/4\pi\varepsilon_0)^2 F_K$, with the King factor
$F_K$ the source used (1 for atoms and spherical tops); the observation is electric and
includes the vibrational part. Measured cross sections determine $\alpha$ to 1-10 %,
against 0.01-1 % from refraction, so they are validators, not fitted (TC-150).

Refraction measures $\alpha + \chi_0$, with $\chi_0$ the magnetic susceptibility per
particle (`[model] chi0_au`); theory and dielectric data give $\alpha$ alone, and are marked
`electric = true` (see the conventions page). For helium, refraction also carries the
retardation term $c\,\omega^2$ (`[model] retardation_au`), added in the same way (TC-705).

A molecule's polarizability depends slightly on temperature. A fit that declares
`[model] reference_temperature_k` brings a dataset measured elsewhere to it with a
`temperature = { measured_k, coefficient, sigma }` option: the values are scaled by
$1 + c\,(T_\mathrm{ref} - T)$, and $\sigma_c\,|T_\mathrm{ref} - T|$ is added to their
relative uncertainty. Validators take the same option (TC-730).

## Reading the fit gallery

The gallery's `fit_overview` plots direct polarizability data against each stored fit
(TC-729). All three rows share a logarithmic vacuum-wavelength axis in nm; static
values are evaluated at zero frequency and displayed at 90 % of the upper wavelength
limit. Positive deviations mean the data or comparison model lies above the fit.

- **Top:** direct datasets, converted to per-particle polarizability. Filled circles
  are fitted data, open squares are validators, and grey dots are other data.
  Configured datasets use the fitter's wavelength windows, additive corrections and
  electric/refractive and infrared conventions; fitted scales are removed. Their
  error bars use the same one-sigma uncertainty as the observation, including shape
  uncertainties and correction allowances. Other datasets retain their stated
  uncertainty and convention. Datasets marked `check` or `excluded` are omitted.
- **Middle:** percentage deviations of published formulas (solid), legacy models
  (dashed) and database comparison models (dotted) from the fitted refractive
  polarizability, each over its own stated range.
- **Bottom:** normalised residuals of the fitted direct observations, using exactly
  their fitting uncertainties. These have the opposite sign to the optimiser's
  residuals, so that positive still means data above fit.

To put the configured data on the same refractive reference as the shaded band, let
$a_i$ and $\sigma_i$ be the corrected observation and its uncertainty, $m_i$ the
model's prediction in that observation's convention, $s_d$ its fitted scale (1 for
validators), and $\alpha_{\rm ref}$ the full refractive polarizability. The top row
shows

$$
100\,\frac{a_i/s_d-m_i}{\alpha_{\rm ref}},
\qquad \text{with error bar}\quad
100\,\frac{\sigma_i/s_d}{|\alpha_{\rm ref}|}.
$$

This is equivalent to adding $\alpha_{\rm ref}-m_i$ to the observation after removing
its fitted scale. For example, an electric datum receives the model's magnetic and
retardation terms; electronic-only theory also receives its infrared contribution.
The bottom row shows $(a_i-s_d m_i)/\sigma_i$.

The first two rows shade the fit's propagated one-sigma covariance band about zero.
The vertical limits are fixed at respectively +/-1 %, +/-2 % and +/-6 sigma; values
outside these limits are clipped. Cauchy moments, phase-matching ratios, Rayleigh
cross sections and priors are not plotted here, so the bottom row is only part of the
fit's residual vector. The title's reduced chi-squared comes from the complete stored
fit.

## Robustness

A leave-one-dataset-out analysis refits the configuration without each dataset,
starting from the full solution, and measures the shift in units of the full fit's
uncertainty. It identifies datasets that determine a model alone: a deliberately
chosen theoretical reference, or a single metrology value that fixes the level
(TC-101). It tests sensitivity to dataset selection, not every source of model error.
See {doc}`accuracy` and the {ref}`robustness gallery <sphx_glr_auto_gallery_accuracy_plot_robustness.py>`.
Commands, stored-file locations and reproduction checks are in
{doc}`/developers/fitting`.

## Compact infrared terms

A molecular fit carries its vibrational bands as fixed poles: the Gaussian rules of every
band system (OEP 3), 4187 poles for CO2 and 6059 for N2O. Such a model cannot be printed
or typed into another code. Its compact model, `compact/<gas>`, keeps the fit's
electronic terms, covariance and rotational system and replaces the vibrational poles by
a few effective ones $(F_j, \Omega_j)$, $F_j \ge 0$, chosen to minimise

$$
\sum_{i}\sum_{m\in\{0,2,3\}}
\left(\frac{D_m[\alpha^{K}_\mathrm{IR}](\omega_i) - D_m[\alpha_\mathrm{IR}](\omega_i)}
{\sigma_m(\omega_i)}\right)^2
\quad\text{subject to}\quad \sum_j F_j = M_0,\;\; \sum_j \frac{F_j}{\Omega_j^2} = M_{-2},
$$

with $D_0[\alpha] = \alpha$ and $D_m[\alpha] = \partial^m(\omega\alpha)/\partial\omega^m$
(proportional to the group-velocity and third-order dispersion in the dilute limit), and
$M_0$ and $M_{-2}$ the vibrational poles' optical-side limit and static contribution, kept
exactly. The weights $\sigma_m$ are the fit's absolute statistical one-sigma of each
quantity, $\sqrt{g_m^\mathsf{T}\Sigma g_m}$ from its covariance with the Birge scaling
(`FittedModel.dispersion_sigma_au`); they stay finite where $D_2$ crosses zero. The
strengths follow exactly for given energies, by an equality-constrained non-negative
least-squares solve over every support; the energies are searched from fixed starting
points (a Gaussian rule, the strongest bands, two spreads and seeded random starts).
The method is OEP 2's ({doc}`oeps/archive/oep-0002-compact-infrared`); the compact models are
derived from the stored fits and never fitted to data.

**Where a compact model holds.** On the fit's qualified intervals inside its optical
range (`range_nm`); the infrared windows beyond are left to the fit, and evaluating the
compact model there warns and names the fit. Every compact pole is bounded to an
absorbing region, below, between or above those intervals, between the region's lowest and
highest parent pole, so none can lie inside an interval. CH4's and N2O's intervals end
short of their optical ranges: the fits' narrow absorbing gaps from 1122 nm (CH4) and
1516 nm (N2O) put band-rule poles within a few cm$^{-1}$ of the interval edges, where no
few poles can follow them (TC-800); the full models hold beyond.

**The tolerance.** Every residual on a check grid four times as fine as the fitted one
(1600 points per interval, geometric in wavelength) is below 0.1 of the fit's sigma in
$\alpha$, $D_2$ and $D_3$; every datum of the fit's datasets and validators in the
compact intervals, or at shorter wavelengths, moves by less than 0.1 of its own sigma;
the static values are the fit's to $10^{-10}$. The tolerance is a dilute-limit one; at
10 bar the relative error in $n - 1$ is the relative error in $\alpha$ to $10^{-3}$ of
itself. The reduction's worst residuals, in the fit's sigma (TC-799):

| Model | Compact interval | Vibrational poles of the fit | Compact poles | $\alpha$ | $D_2$ | $D_3$ |
|---|---|---:|---:|---:|---:|---:|
| `compact/CO` | 167-2000 nm | 96 | 3 | 0.0017 | 0.0057 | 0.049 |
| `compact/CO2` | 250-1400 nm | 4187 | 5 | 0.0022 | 0.0015 | 0.046 |
| `compact/CH4` | 240-1070 nm | 943 | 5 | 0.0007 | 0.00059 | 0.01 |
| `compact/N2O` | 320-1500 nm | 6059 | 5 | 0.00056 | 0.00021 | 0.069 |
| `compact/SF6` | 230-5000 nm | 3 | 3 | 0 | 0 | 0 |

SF6's three cross-section poles are already compact and are carried through unchanged.
One pole fewer fails in every other case, worst in $D_3$. The compact poles are effective
ones, not band centres.

```python
from optmat.gas.models.registry import default_model, model

co2 = model("compact/CO2")
assert co2.parent == "fit/CO2"
assert default_model("CO2").name == "fit/CO2"  # never a default
assert max(co2.residuals_sigma.values()) < 0.1
```

**Regenerating.** See {doc}`/developers/fitting` for configuration fields, commands,
checks and the required order after a refit.

**For print.** `optmat.gas.models.compact.coefficients(gas)` returns the polarizability
form, the terms' $f_k$ and $\lambda_k$ in descending energy with the qualified intervals
and the electronic covariance in those parameters; with `form="sellmeier"` and a state,
the dilute-limit $B_k$ and $C_k$ at that state, with their error against the full model.
`format_coefficients` prints either to eleven significant digits, which keep the rounding
below 0.01 of the fit's sigma (TC-801). {doc}`gas_models` has both tables for each
molecule.

The compact models are what [OEP 2](oeps/archive/oep-0002-compact-infrared.md)
proposed and implemented: its theory, decisions and phase reports are the record.
