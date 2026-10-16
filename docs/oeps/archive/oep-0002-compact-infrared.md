# OEP 2 — Compact infrared terms for the molecular gases

| | |
|---|---|
| **OEP** | 2 |
| **Title** | Compact infrared terms for the molecular gases |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | GPT-6 (AI); John Travers (human) |
| **Type** | Feature |
| **Created** | 2026-09-29 |
| **Resolution** | 2026-09-30: accepted by the author, John Travers, after the human and AI reviews; merged directly to `main` |
| **Requires** | [OEP 3](oep-0003-infrared-windows.md): the parents' qualified intervals and representation bound |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phases 1 and 2 merged from `feature/oep-0002-compact` |
| **Manual pages** | {doc}`../../fitting` (*Compact infrared terms*), {doc}`../../gas_models` (*The compact coefficients*) |

## Part I — Proposal

### 1. Abstract

The fitted models of the infrared-active molecules carry their vibrational bands from
HITRAN as fixed poles, one per band: 4296 for CO2, 2518 for N2O, 295 for CH4 and 30 for
CO. The model is correct, but it cannot be printed, typed into another code or checked
by hand. This OEP proposes a **compact model** beside each molecular fit, `compact/<gas>`.
It keeps the fit's electronic terms and covariance unchanged, and replaces the fixed
poles by a few effective infrared poles. These are derived from the stored fit by a
constrained least-squares fit to its infrared contribution over every interval in which
the fit is qualified. The fit keeps two moments of the bands exactly: the optical-side
limit and the static contribution.

This OEP is sequenced after OEP 3, which extends the parent fits into the infrared
windows between their bands, with several qualified intervals per model and a more
accurate representation of each band. The compact models follow whatever intervals the
parents qualify. Each compact pole is kept inside an absorbing region, outside every
qualified interval.

A feasibility calculation against today's parents, each qualified on one interval from
the ultraviolet to 1.4-2 µm, finds the following:
- 2 (CO), 3 (CH4), 4 (N2O) and 6 (CO2) poles reproduce the full model everywhere in its
  range to within 0.06 of its own statistical uncertainty;
- the check covers the polarizability, the group-velocity dispersion and the third-order
  dispersion;
- each molecule then needs 4 to 10 Sellmeier-type terms in all;
- the reproduction extends to every dataset and validator of the parent fits inside that
  range, and to its static values, which the constrained moment preserves exactly.

With OEP 3's windows the number of poles will grow with the band systems that bound
them; Phase 1 recomputes it.

The compact models are derived and never fitted to data. They are never a default. They
are regenerated and tested whenever their parent fit changes. Success is judged by one
stated tolerance, met at every wavelength of every qualified interval. The proposal also
defines how the coefficients are exported for print:
the polarizability poles with the density conversion, and, as a convenience, an $n - 1$
Sellmeier form at a stated state, with its error against the full model stated.

### 2. Motivation and scope

**The use.** Three readers need the molecules' dispersion in a closed form:

- **a publication** of optmat's gas models, whose tables must let a reader compute
  $n(\lambda)$ without the package;
- **propagation codes** that do not call optmat at run time and embed a formula;
- **anyone checking a number by hand** against a published formula.

The alternatives in the literature for these gases have no infrared term. He et al.'s
UV-visible formulas for CH4, SF6 and N2O [10] overstate SF6's group-velocity dispersion at
800 nm by 55 % and miss Shelton's near-infrared phase-matching ratios by 3-28 %
(`docs/gas_models.md`, TC-104). Air formulas with infrared terms exist, but not for the
pure molecular gases (§3).

**In scope:**
- compact models for CO, CO2, CH4 and N2O, derived from their stored fits and qualified
  over each fit's qualified intervals (`Model.intervals`, OEP 3; today a single interval,
  `Model.wavelength_range`);
- SF6, whose three poles are already compact, carried through unchanged so that every
  molecule has the same interface;
- the derivation tool, the stored compact models and their tests;
- an export of the coefficients in two printed forms, and a manual section with the
  tables.

**Out of scope:**
- **Qualifying the parents in the infrared windows.** That is OEP 3.
  - Today's parents declare one interval each; CO2's is 250-1400 nm. Its 8-12 µm window
    appears only in its configuration's description and on the accuracy page, and is
    not qualified.
  - There, the parent's one pole per band puts 224 poles among the lines, which alone
    contribute a median of 2.5 times its statistical sigma (§8). No few-pole model can
    follow that, and the parent is not accurate there either.
  - OEP 3 replaces those poles, folds the weak in-window bands and qualifies the
    windows. The compact models then cover them.
  - Until then, the compact CO2 model warns in the window, like the parent, and differs
    from it by up to 3 % in the polarizability.
- **Data wavelengths outside every qualified interval.** Some parent fits use data
  beyond their declared range today: CO2 uses Old et al.'s points to 1529 nm and Egan
  et al.'s 1542 nm, and CH4 uses Zhao's 1566 nm. The compact model reproduces these as
  reported numbers (§8), without a promise. After OEP 3 most of them lie in qualified
  windows.
- **Ozone.** Its model takes the measured absorption by Kramers-Kronig and is complex
  (`docs/gas_models.md`); it is not a sum of fixed poles.
- **Changing the full models.** The compact model represents the stored fit. It does not
  refit the electronic terms or the HITRAN strengths. The band representation is OEP 3's
  (§11, §13 Q3).
- **The absorption region.** Like the full model, the compact one describes the index
  away from the bands. It is not a model of the infrared absorption.
- **Temperature dependence.** HITRAN's strengths refer to 296 K, and the molecular
  polarizabilities to 293.15 K (TC-730), as in the full models.
- **Air and water vapour.** optmat has no air or water model (OEP 1 proposes an air
  model; the compact terms would compose with it).

### 3. Context

**What optmat does today.**
- {mod}`optmat.gas.infrared` reads the band totals of
  `data/external/hitran/<gas>_bands.csv`. Per band these are the sum of the line
  intensities $S$ and of $S/\tilde\nu^2$, and the span of its lines (`min_cm`, `max_cm`).
  Each band becomes one pole that keeps the band's total strength and its static
  contribution. Pure rotational bands are left out (TC-156). See `docs/data.md`,
  *Infrared bands (HITRAN)*.
- A fit configuration's `[model] infrared_max_wavenumber` adds the bands up to that
  wavenumber as fixed trailing terms of the fitted model, outside its covariance
  (`docs/fitting.md`, {class}`~optmat.gas.models.fitted.FittedModel`, `n_fixed`).
- The molecular fits (`docs/gas_models.md`), as stored on 2026-09-30:

  | Fit | Electronic terms | Fixed poles | Highest pole | Declared range | Its lowest frequency |
  |---|---|---|---|---|---|
  | `fit/CO` | 2 | 30 | 4262.0 cm⁻¹ | 167-2000 nm | 5000 cm⁻¹ |
  | `fit/CO2` | 4 (2 lines, 2 effective) | 4296 | 6996.7 cm⁻¹ | 250-1400 nm | 7142.9 cm⁻¹ |
  | `fit/CH4` | 2 | 295 | 4499.7 cm⁻¹ | 240-1400 nm | 7142.9 cm⁻¹ |
  | `fit/N2O` | 2 | 2518 | 4499.5 cm⁻¹ | 320-2000 nm | 5000 cm⁻¹ |
  | `fit/SF6` | 2 | 3 (cross-sections) | 1577.5 cm⁻¹ | 230-5000 nm | 2000 cm⁻¹ |

  Every declared range lies above every fixed pole. The closest case is CO2: its highest
  pole is 2.0 % below the range's lowest frequency. OEP 3 changes this inventory: every
  band of the tables kept, several poles per band, and several qualified intervals with
  bands between them.

**Register entries.**
- TC-104: the HITRAN strengths are checked against Shelton's phase-matching ratios.
- TC-146: the CH4, SF6 and N2O fits with the infrared bands.
- TC-148: the CO2 fit and its range.
- TC-156: N2O's rotational (orientation) term is left out.
- TC-720: Voronin & Zheltikov's generalized Sellmeier equation for air, the closest
  precedent (below).
- TC-721 and TC-722: Mathar's line-by-line models of humid air.
- TC-723: Tang's bounds on one-term formulas, and what they say of the per-band pole.
- TC-725 and TC-726: the refits of `fit/CO2`, and its electronic form (two weak bands as
  lines, two effective oscillators from 10 eV).
- TC-730: the molecules' reference temperature, 293.15 K.

**Prior work** (curated as method records):
- **Voronin & Zheltikov** [1] represent air's infrared (HITRAN 2012, about 650,000 lines)
  by two poles per band, 15 bands and 30 poles in all. Only each band's total strength is
  constrained; the pole wavelengths are fitted. $n - 1$ is linear in the partial
  densities, without Lorentz-Lorenz, over 0.3-13 µm. Their CO2 terms are the only
  separate infrared terms for CO2 in the literature. Against optmat's HITRAN 2024 bands
  they are 1.7 % low in strength and 1.8 % low in static contribution, a match to the
  principal isotopologue alone (TC-720; an inference, not stated by them). Their
  per-species electronic terms are tuned as a set to Ciddor's air formula and are not
  per-species polarizabilities (argon 6.5 % low; TC-720).
- **Mathar** sums every HITRAN line of water and CO2, with ultraviolet pseudospectra and
  no fitting [2]. His later fits of humid air are polynomials in wavenumber, one per
  window [3]. He reports that a two-resonance Sellmeier form fitted better at a fixed
  state, but that expanding its parameters in temperature, pressure and humidity "did
  not achieve competitive stability" ([3], p. 473; TC-722).
- **Tang** [4] proves bounds on one-term dispersion formulas for positive oscillator
  distributions; Tang & Poon [5] generalise them to several moments. §4 applies them to
  the infrared poles.

**What is new here.**
- Per molecule rather than per mixture.
- A polarizability with optmat's density conversion (Lorentz-Lorenz, real-gas density,
  refractivity virial), rather than $n - 1$ linear in density.
- Few poles for all the bands together, fitted on the qualified intervals, with two moments
  exact and the error stated against the model's own uncertainty.
- All isotopologues, and electronic terms fitted to measurements with a covariance.

Mathar's instability does not arise: the poles are properties of the molecule, and the
state enters only through the density.

**Related OEPs.** OEP 3 (Draft) extends the parents into the infrared windows and is to
be implemented first (§12). OEP 1 (Draft) declares the propagation-code interface.

**Roadmap.** None of the informal items covers this; the draft is listed under *Draft
enhancement proposals*.

### 4. Theory

**The infrared terms.** In SI, the electronic and vibrational oscillators of a molecule
give the dynamic dipole polarizability

$$
\alpha(\omega) = \frac{e^2}{m_e}\sum_k \frac{f_k}{\omega_k^2 - \omega^2}
$$

(C m² V⁻¹; in atomic units $\alpha = \sum_k f_k/(E_k^2 - \omega^2)$, the form the code
uses), with dimensionless oscillator strengths $f_k$. A vibrational line of integrated
intensity $S$ (HITRAN's definition, cm/molecule at 296 K, which includes the lower-state
population and stimulated emission) has $f = S/(\pi r_e)$, with $r_e$ the classical
electron radius in cm (`optmat.gas.infrared`). Write the infrared part as
$\alpha_\mathrm{IR}(\omega) = (e^2/m_e)\sum_{b} f_b/(\omega_b^2 - \omega^2)$, summed over
the $B$ fixed poles of the stored fit, and its moments

$$
M_k = \sum_b f_b\,\omega_b^{k}.
$$

Above every band ($\omega > \max_b \omega_b$, which holds throughout today's declared
ranges, §3)

$$
\alpha_\mathrm{IR}(\omega) = -\frac{e^2}{m_e}\left(\frac{M_0}{\omega^2}
+ \frac{M_2}{\omega^4} + \frac{M_4}{\omega^6} + \cdots\right),
$$

and below every band $\alpha_\mathrm{IR}(\omega) = (e^2/m_e)(M_{-2} + M_{-4}\omega^2 +
\cdots)$. $M_0$ is the vibrational Thomas-Reiche-Kuhn sum and fixes the optical-side limit.
$(e^2/m_e) M_{-2}$ is the vibrational static polarizability. In OEP 3's windows, bands lie
on both sides of the frequency, and neither expansion holds. The reduction below does
not rely on either: it fits the parent's terms directly on the qualified intervals.

**The parent's representation of a band.** Today one pole per band, a lower bound.
OEP 3 replaces it by a Gaussian rule of a few poles per band, built from the lines. The
reduction below represents whatever the parent stores, and inherits its representation
without adding to it. Today's pole for band $b$ takes
$f_b = \sum_{\ell\in b} f_\ell$ and $\omega_b^2 = \sum f_\ell / \sum (f_\ell/\omega_\ell^2)$,
keeping the band's $M_0$ and $M_{-2}$. That is Tang's construction [4] applied to one
band: below the band the pole bounds the band's contribution from below, and above it
the pole bounds its magnitude from below. Above the band, the pole's magnitude falls
short of the band's by, to leading order,

$$
\frac{e^2}{m_e}\,\frac{M_2^{(b)} - \left(M_0^{(b)}\right)^2/M_{-2}^{(b)}}{\omega^4} \;\ge\; 0 ,
$$

non-negative by the Cauchy-Schwarz inequality. For CO2's 4.3 µm band
it is about $6\times10^{-6}$ of the band's contribution at 633 nm (TC-723). At the edges
of the infrared windows the error is far larger, which is why OEP 3 replaces it.

**The compact representation.** Replace the $B$ poles by $K \ll B$ poles
$(F_j, \Omega_j)$, $F_j \ge 0$, chosen to minimise

$$
\chi^2_\mathrm{red} = \sum_{i}\sum_{m\in\{0,2,3\}}
\left(\frac{D_m[\alpha_\mathrm{IR}^{K}](\omega_i) - D_m[\alpha_\mathrm{IR}](\omega_i)}
{\sigma_m(\omega_i)}\right)^2
\quad\text{subject to}\quad \sum_j F_j = M_0,\;\; \sum_j \frac{F_j}{\Omega_j^2} = M_{-2},
$$

over a grid $\omega_i$ spanning every qualified interval, where $D_0[\alpha] = \alpha$ and
$D_m[\alpha] = \partial^m(\omega\alpha)/\partial\omega^m$ for $m = 2, 3$. These are
proportional to the group-velocity dispersion $\beta_2$ and the third-order dispersion
$\beta_3$ in the dilute limit ({meth}`~optmat.gas.models.fitted.FittedModel.dispersion_sigma_relative`).
The two equalities keep the optical-side limit and the static contribution exactly, so
the compact model has the parent's static polarizability. The strengths enter linearly.
For given $\Omega_j$ they follow from an equality-constrained non-negative linear
least-squares problem (variable projection, §7). Only the $K$ energies are searched,
each within one absorbing region (below).

**The weights.** $\sigma_m(\omega)$ is the parent's absolute one-sigma statistical
uncertainty of $D_m$ of the full polarizability, $\sqrt{g_m^\mathsf{T}\,\Sigma\,g_m}$,
with $g_m$ the gradient of $D_m$ with respect to the fitted parameters and $\Sigma$ their
covariance. The fixed infrared terms carry no statistical uncertainty.
- OEP 3's representation bound on the parent's infrared terms
  (`representation_bound_au`) is a deterministic bound, not a standard deviation. It is
  not combined with $\sigma_m$.
- The compact model reproduces its parent. Its error against OEP 3's line-sum reference
  is at most the parent's bound plus the compact model's own residual, and both are
  stated.
- OEP 3's rotational band system is carried unchanged, outside the reduction. The
  moments $M_0$ and $M_{-2}$ constrained here are those of the vibrational systems. The covariance is scaled by the reduced chi-squared where that
exceeds one, the Birge convention of the parent's own uncertainties (`birge=True`, the
default of {meth}`~optmat.gas.models.fitted.FittedModel.polarizability_sigma_au`). The
absolute sigma is finite where $D_2$ crosses zero: at SF6's zero-dispersion wavelength,
1271.2 nm, it is 0.0284 au, against a largest $|D_2|$ of 2576 au over the range. No
floor is needed. For $m \ge 2$ this sigma equals the package's relative uncertainty times
$|D_m|$ wherever $D_m \ne 0$.

**Where the compact model is valid.** The compact model is qualified on the parent's
qualified intervals. The frequencies outside them divide into absorbing regions:
- below the lowest interval;
- between two adjacent intervals, where a band system absorbs.

Only regions that contain parent poles count. Each $\Omega_j$ is bounded to lie within
one such region, between its lowest and highest parent pole, and $K$ is divided among
the regions. So no compact pole can fall inside a qualified interval. This follows from
the bounds, independent of the sampling grid, and the tests assert it directly.

On today's parents there is one region, below the single interval: every declared range
lies above $\max_b \omega_b$ (§3).

**The density conversion.** The tolerance is stated on the polarizability and its
derivatives, which do not depend on the density; $\beta_2$ and $\beta_3$ are
proportional to $D_2$ and $D_3$ in the dilute limit. At a finite density the parent and
the compact model pass through the same conversion (Lorentz-Lorenz, the real-gas
density, the refractivity virial). To first order in the difference, the relative
error in $n - 1$ is then the relative error in $\alpha$, and the ratio of error to
sigma is unchanged, up to a relative correction of the order of the susceptibility.
§8 measures it: $7\times10^{-5}$ at 1 bar and $8\times10^{-4}$ at 10 bar, of itself.
The tolerance is therefore labelled a dilute-limit tolerance, and a test repeats the
comparison at 10 bar.

**Why not moment matching alone.** A $K$-point Gaussian quadrature of the band
distribution in the variable $x = \omega_b^2$, with weights $f_b$, has positive weights
and matches $2K$ moments in $x$: $M_0, M_2, \ldots, M_{4K-2}$, the coefficients of the
optical-side expansion above. Its error is governed by the next term of that expansion
in $(\omega_b/\omega)^2$. The expansion converges slowly at the long-wavelength end of a
range: CO2's highest pole, at 6996.7 cm⁻¹, lies 2.0 % below the range's end at 1.4 µm
(7142.9 cm⁻¹). The feasibility calculation (§8) confirms this. A 6-node rule on the
current CO2 parent misses by 23 sigma in $D_2$ and 7320 sigma in $D_3$, and a 10-node
rule by 1.1 and 469 sigma. Moment matching is kept as a starting point of the search.

**The printed forms.** The compact polarizability is a finite sum of poles, so the
package's density conversion applies unchanged ({doc}`../../conventions`). For print, a
dilute-limit Sellmeier form at a stated state $(p, T)$ follows from
$f/(\omega_k^2 - \omega^2) = (f/\omega_k^2)\,\lambda^2/(\lambda^2 - \lambda_k^2)$:

$$
n - 1 \approx \frac{N\alpha}{2\varepsilon_0}
= \sum_k B_k \frac{\lambda^2}{\lambda^2 - C_k},\qquad
B_k = \frac{N e^2 f_k}{2\varepsilon_0 m_e \omega_k^2},\quad C_k = \lambda_k^2 ,
$$

with $\lambda$ the vacuum wavelength and $N$ the real-gas number density at the stated
state ({mod}`optmat.gas.eos`). The infrared terms have $\lambda < \lambda_k$ and
contribute negatively. The form drops Lorentz-Lorenz, whose leading relative correction
to $n - 1$ is $(n-1)/6$, and the refractivity virial. The export does not rely on that
estimate. It states the form's error as the largest relative difference in $n - 1$ and
in $\beta_2$ from the full parent model at that state, over the qualified intervals,
computed directly. For the four molecules at 101325 Pa and 293.15 K the difference in $n - 1$ is
0.5 to $1.2\times10^{-4}$, the same order as the parent's statistical sigma (§8). The
Sellmeier form is therefore a convenience at its state. The polarizability form carries
the stated tolerance.

**Conventions.** SI throughout, vacuum wavelengths, and the refractive polarizability as
defined on the conventions page: electric plus $\chi_0$, with $\chi_0 = 0$ for the
molecules here. The effective infrared poles are not band centres, and are not to be
read as spectroscopy.

### 5. Primary sources and data

- **No new data.** The compact models derive from the stored fits
  (`data/gas/models/fitted/<gas>.json`), which already contain the HITRAN 2024 [6] band
  poles (`data/external/hitran/`, TC-104, TC-148). No transcription, extraction or
  vocabulary quantity is needed.
- **The tolerance is set by the parent fit's covariance** (§4, §7), not by external
  data.
- **Reproduction checks.** The compact model is a representation of its parent, so its
  checks test the reproduction, not the parent's physics.
  - Every fitted dataset and validator of the parent is re-evaluated with the compact
    infrared terms, through the fit's own observation code, at its own wavelengths
    (harmonics included). The parent's fitted scales and partner models stay fixed.
  - The static validators test the constrained $M_{-2}$: the CH4 and CO2 dielectric
    values and CO's total static polarizability.
  - Marchetti & Simili's CO2-laser values [7] and Shelton's ratios [8, 9] enter the
    parent fits, so they are reproduction checks, not independent validation.
- **Validation of the physics** is the parent's (its validators, TC-104, TC-146,
  TC-148, and OEP 3's infrared validators). The compact model inherits it within the
  qualified intervals, and §8 shows the difference there is negligible against every
  datum's uncertainty.
- **Methods cited:** Voronin & Zheltikov [1], Mathar [2, 3], Tang [4] and Tang & Poon [5],
  all curated.
- **Uncertainty.** The compact model carries the parent's covariance of the electronic
  terms. Its infrared poles are fixed, as the HITRAN poles are. The reduction error is
  stated per gas as the largest residual, in units of the parent's sigma, of $\alpha$,
  $D_2$ and $D_3$ over the qualified intervals. It is recorded in the stored compact model
  and quoted in the manual.

### 6. Proposed interface

- **Models.** `compact/<gas>` for CO, CO2, CH4, N2O and SF6.
  - Each is loaded by `optmat.gas.models.compact.compact_models()` and joins
    `optmat.gas.models.registry` beside the fitted, published and comparison models.
    It is listed by `registry.models(gas)` and returned by
    `registry.model("compact/CO2")`, and never by `registry.default_model`.
  - It is a {class}`~optmat.gas.models.fitted.FittedModel` with the parent's leading
    terms and covariance and $K$ fixed trailing terms. So `polarizability`,
    `dispersion_au`, `polarizability_sigma_au` and `dispersion_sigma_relative` work
    unchanged.
  - Its qualified intervals are the parent's: `intervals` (OEP 3) and the envelope
    `wavelength_range`. Outside every interval, `ModelRangeWarning`, as for the parent.
    Today each parent has a single interval.
  - Its record names the parent (`fit/<gas>`), the SHA-256 of the stored parent record,
    the parent's `config_sha256`, the digest of the reduction configuration, the
    tolerance and the achieved residuals.
- **Storage and configuration.**
  - The stored models are `data/gas/models/compact/<gas>.json`, written by
    `tools/fit/compact.py <gas>` from the stored parent.
  - The reduction's settings are a separate file, `data/gas/compact/<gas>.toml`: $K$,
    the constraints, the grid and the tolerance.
  - The parent's configuration is not touched, so its digest, stored fit and robustness
    record do not change (§13 Q11).
- **Export.**

  ```text
  from optmat.gas.models.compact import coefficients
  table = coefficients("CO2")                  # the polarizability form
  table = coefficients("CO2", form="sellmeier", pressure=101325.0, temperature=293.15)
  ```

  - `coefficients` returns a frozen `CompactCoefficients`:
    - the model and parent names and their digests;
    - `form`, and for the Sellmeier form the state (Pa, K) and the number density used
      (m⁻³, from `optmat.gas.eos`);
    - the term `labels`, and `kinds` (line, effective or infrared);
    - for the polarizability form, `strengths` $f_k$ (dimensionless) and `wavelengths`
      $\lambda_k$ (m), and the constant term $\chi_0$ (C m² V⁻¹; zero for these
      molecules);
    - for the Sellmeier form, `B` (dimensionless) and `C` (m²);
    - the `intervals` (m) in which the coefficients are qualified. A printed table
      states them, and says that between them the gas absorbs;
    - the reduction's achieved residuals, and for the Sellmeier form its stated error
      (§4): the largest relative differences in $n - 1$ and in $\beta_2$ from the full
      parent at that state, over those intervals;
    - the covariance of the electronic terms, transformed to the printed parameters.
  - Terms are ordered by descending energy (ascending $\lambda_k$): electronic, then
    infrared.
  - Pressure is in Pa and temperature in K, both within the range where the gas's
    equation of state is valid.
  - A formatter writes the printed tables with 10 significant digits, the number the
    round-trip test requires (§8).
- **Refused.** `compact/<gas>` for a gas without fixed infrared terms (atoms, H2, D2, N2,
  O2) and for O3: `OptmatError`.

### 7. Implementation strategies

**Where.** In `optmat.gas`: a module `optmat.gas.models.compact` (loading, export) and the
tool `tools/fit/compact.py` (derivation). The registry's `_all_models` gains the compact
models. There are no core changes, no vocabulary changes and no new dependencies (NumPy
and SciPy's `least_squares`, already used by the fits).

**The derivation. Recommended: constrained least squares on the qualified intervals.**
- **The inner problem.** For fixed energies $\Omega$, minimise $\lVert A F - b\rVert^2$
  subject to $C F = d$ (the $M_0$ and $M_{-2}$ rows) and $F \ge 0$.
  - Solve it exactly by enumerating the supports $S \subseteq \{1..K\}$. On each support,
    solve the equality-constrained least squares by the null-space method (a particular
    solution of $C_S F_S = d$, plus the least-squares solution in the null space of
    $C_S$).
  - Keep the supports whose solution is non-negative and meets the equalities to a
    relative $10^{-10}$, and take the one of lowest cost.
  - The optimum of this convex problem is the equality-constrained minimiser on its own
    support, so the enumeration finds it. The $2^K$ supports ($\le 1024$ for
    $K \le 10$) cost milliseconds.
  - For larger $K$, which OEP 3's windows may need, an active-set solver of the same
    problem replaces the enumeration. It is held to the same moment tolerance, and
    checked against the enumeration at $K \le 10$.
  - Strengths are scaled by $M_0$ for conditioning.
  - Penalty rows are not used: they do not make a moment exact.
- **Feasibility.** With both equalities, a solution with $F \ge 0$ exists exactly when
  $\min_j \Omega_j \le \bar\omega \le \max_j \Omega_j$, with
  $\bar\omega = (M_0/M_{-2})^{1/2}$. The search guarantees this by its bounds.
  - At least one energy is bounded to parent poles below $\bar\omega$, and one to
    parent poles above it, each within its absorbing region.
  - The others lie anywhere within their regions.
  - For $K = 1$ the pole is fixed at $\bar\omega$ (Tang's pole for all the bands
    together). This is possible only when one region holds every parent pole.
- **The outer problem.** Bounded nonlinear least squares in $\log\Omega$, on the
  weighted residual vector of the inner solution.
- **Starting points** (all tried, with a recorded seed), within each absorbing region
  for its share of $K$:
  - the Gaussian rule in $x = \omega_b^2$ of the region's poles (§4);
  - the energies of the region's strongest bands at least 5 % apart;
  - two log-spaced spreads over the region;
  - for CO2, twelve log-uniform random starts. Without them the $K = 6$ search stopped
    in a local minimum with a pole of zero strength.
- **Selection.** Among the converged candidates, those meeting every tolerance on the
  check grid are preferred, and of those the lowest cost is kept. The lowest cost alone
  can have a worse largest residual and fail a threshold that another candidate passes.
  If none passes, $K$ is raised.
- $K$, and its division among the regions, is the smallest that meets the tolerance.
  It is found region by region and recorded in the reduction configuration.
- The derivation is deterministic: a fixed grid, fixed starts and seed, and a fixed
  selection rule.

**The tolerance. Recommended:**
- the reduction: every residual below 0.1 of the parent's sigma in $\alpha$, $D_2$ and
  $D_3$ (§4), at every point of the check grid over every qualified interval. The sigma
  is the parent's statistical one (§4). In OEP 3's windows the parent's representation
  bound is stated beside it;
- rounding: the printed coefficients, read back from their decimal strings, change the
  same quantities by less than 0.01 sigma;
- reproduction: every datum of the parent's datasets and validators whose wavelengths
  all lie in qualified intervals, or at shorter wavelengths than the first (harmonics in
  the ultraviolet), moves by less than 0.1 of its own sigma. Its chi-squared change is reported rather than bounded, since that change,
  $2 r^\mathsf{T}\delta + \delta^\mathsf{T}\delta$, depends on the parent's residuals
  $r$;
- static values: equal to the parent's to a relative $10^{-10}$ (the constrained
  $M_{-2}$).

**Alternative:** a fraction of the working estimate. It is looser and not computable,
since the working estimates are judgements.

**The grid.** Geometric in wavelength over each qualified interval: 400 points fitted per
interval, and the residual checked on four times as many. The per-datum checks of §5 add the data wavelengths.

**Option B, moment matching alone.** Deterministic and parameter-free, with guaranteed
positivity. It fails the tolerance near the long-wavelength end (§4, §8), and has no
meaning in a window with bands on both sides. It is kept as a starting point, not the
method.

**Option C, a joint refit.** Fit the compact infrared poles together with the electronic
terms to the data. Rejected: it makes a second, different model with its own
uncertainty. The compact model should be a faithful representation of the published fit.

**Cost.** On today's parents, seconds per gas for CO, CH4, N2O and SF6. CO2, with its
4296 poles on the grid and its extra starts, takes about 30 s. OEP 3's parents have
more poles and more intervals: minutes per gas are acceptable, since the derivation runs
only when a parent changes. Checking a stored compact model's residuals and
its parent's digest takes well under a second and runs in the default test lane. The
re-derivation runs in the `fits` lane.

### 8. Verification and validation

- **Feasibility on today's parents** (computed 2026-09-30 against the stored fits after
  TC-730, with the method and tolerance of §7; one qualified interval each). It shows
  that the method works and what it costs. Phase 1 recomputes it on OEP 3's parents,
  with their windows. The last three columns are the worst residual over the declared
  range, on the check grid, in units of the parent's sigma:

  | Gas | Parent terms | $K$ | Poles (cm⁻¹) | Terms in all | $\alpha$ | $D_2$ | $D_3$ |
  |---|---|---|---|---|---|---|---|
  | CO | 2 + 30 | 2 | 2144.29, 4261.20 | 4 | 0.0005 | 0.006 | 0.057 |
  | CH4 | 2 + 295 | 3 | 1309.87, 3026.07, 4290.76 | 5 | < 0.0001 | < 0.0001 | 0.0002 |
  | N2O | 2 + 2518 | 4 | 999.58, 2222.49, 3521.16, 4410.52 | 6 | 0.0005 | 0.0002 | 0.006 |
  | CO2 | 4 + 4296 | 6 | 718.72, 2383.86, 4140.41, 6288.96, 6916.97, 6968.31 | 10 | 0.0002 | 0.0002 | 0.011 |
  | SF6 | 2 + 3 | 3 | 615.60, 945.05, 1577.53 (the parent's, unchanged) | 5 | 0 | 0 | 0 |

  - One pole fewer fails in every case, worst in $D_3$: CO $K = 1$ by 2700 sigma, CH4
    $K = 2$ by 2.2, N2O $K = 3$ by 2.6 and CO2 $K = 5$ by 0.26.
  - $M_0$ holds to a relative $10^{-11}$ and $M_{-2}$ to $2\times10^{-15}$. No energy
    lies on a bound, and every energy lies below its range's lowest frequency (§3).
  - Parents (SHA-256 of `data/gas/models/fitted/<gas>.json`):
    - CO `bebaa3c0a2524636469a306ca6134f3f04813887796b8662fdc9fc9040b2f9be`;
    - CH4 `45136078db914ef910dd5d60ff2e3e762ebb7e27f38f6ff2276c4e848193ebbd`;
    - N2O `89b74841e95bbf7fa9e7c4457fdbf15a6a1bc8a1d0a06e410de5ac35ae914576`;
    - CO2 `aef7f296e07513017a19d99e04f1a0cf95053105026dc60d4da17c14b881f17b`;
    - SF6 `91a123f67c5a599d499672e58508db7235e11963606a9404ee6d21773d5b88c5`.
  - The infrared coefficients, $F_j$ and $\Omega_j$ in hartree, to evaluate the residuals
    independently:

    | Gas | $F_j$ | $\Omega_j$ (hartree) |
    |---|---|---|
    | CO | 1.165240679e-05, 8.649074412e-08 | 9.770118942e-03, 1.941543689e-02 |
    | CH4 | 6.018792862e-06, 1.296480419e-05, 8.835528366e-07 | 5.968208577e-03, 1.378776856e-02, 1.955015025e-02 |
    | N2O | 1.051371077e-05, 6.867801594e-05, 2.361823676e-06, 7.821345467e-08 | 4.554409104e-03, 1.012639765e-02, 1.604360791e-02, 2.009579522e-02 |
    | CO2 | 1.232725699e-05, 1.141932460e-04, 1.340842468e-06, 1.597703545e-09, 1.450522007e-10, 1.782683831e-09 | 3.274751732e-03, 1.086165586e-02, 1.886509979e-02, 2.865459354e-02, 3.151604114e-02, 3.174995254e-02 |

  - The sigmas were checked against the package's own
    `polarizability_sigma_au` and `dispersion_sigma_relative` times $|D_m|$, which they
    equal.
- **Reproduction of the parents' data** (the compact infrared terms in the fit's own
  observation code, fitted scales and partner models fixed):
  - inside each declared range, the largest change of any datum is 0.0002 of its sigma
    (CO 0.00015, CH4 0.00002, N2O 0.00018, CO2 0.00017);
  - the static validators are equal to a relative $4\times10^{-16}$;
  - beyond the range (reported, not asserted):
    - CH4 at 1566 nm, 0.00003 sigma;
    - CO2 at 1476-1529 nm (Old et al.), 1.3 sigma, and at 1542 nm (Egan et al. 2026),
      0.15 sigma;
    - CO2 at 9.2-10.7 µm (Marchetti, Willetts, Simmons), 1.2 to 10 sigma, up to 3 % in
      $\alpha$.
- **CO2's window on today's parent** (the grounds of §2, and one motivation of OEP 3):
  224 parent poles lie in 833-1250 cm⁻¹, carrying
  $1.45\times10^{-5}$ of $M_0$. The lines of 258 kept bands reach into the window and
  together cover all of it. On a 4001-point grid, those 224 poles alone contribute a
  median of 2.5 sigma, and more than 0.1 sigma at 99.5 % of the points.
- **Rounding.** Coefficients rounded to $d$ significant digits (the electronic and
  infrared $f_k$ and $\lambda_k$) change the worst residual by the following amounts.
  - CO2: 16 sigma at 6 digits, 0.32 at 8, 0.034 at 9 and 0.0034 at 10. $D_3$ near
    1.4 µm needs the digits.
  - CO, CH4 and N2O: below 0.01 sigma from 7 or 8 digits.
  - So the printed tables use 10 significant digits.
- **Finite density.** At 293.15 K and 1 and 10 bar, the relative error in $n - 1$ of the
  compact model against the parent equals the relative error in $\alpha$ to
  $7\times10^{-5}$ and $8\times10^{-4}$ of itself (§4).
- **The Sellmeier form** at 101325 Pa and 293.15 K, against the full parent (with
  Lorentz-Lorenz, the real-gas density and the virial). The relative error in $n - 1$
  is:

  | Gas | Relative error in $n - 1$ |
  |---|---|
  | CO | $-8.7$ to $-5.1\times10^{-5}$ |
  | CH4 | $-1.2$ to $-1.1\times10^{-4}$ |
  | N2O | $-8.4$ to $-7.6\times10^{-5}$ |
  | CO2 | $-8.8$ to $-8.0\times10^{-5}$ |

  This is mainly Lorentz-Lorenz. At 273.15 K it is about 8 % larger.
- **Invariants:**
  - $M_0$ and $M_{-2}$ are exact to a relative $10^{-10}$;
  - every $F_j \ge 0$;
  - every $\Omega_j$ lies within its absorbing region, outside every qualified interval:
    an analytic check, independent of the grid;
  - the compact model equals the parent at every electronic term and in its covariance.
- **Freshness:** a default-lane test checks each stored compact model's parent SHA-256
  and reduction digest against the files, and its recorded residuals on the check grid.
  A `fits`-lane test re-derives each compact model and compares it with the stored one,
  as `test-fits` does for fits.
- **Tolerance:** a test evaluates the residuals on the check grid and asserts the
  tolerance of §7, and repeats the $n - 1$ comparison at 10 bar.
- **Reproduction:** a test asserts the per-datum and static tolerances of §7, through
  the fit's observation code at the data's own wavelengths, harmonics included.
- **Round trip:** a test formats both printed forms, parses the decimal strings back and
  asserts the rounding tolerance. For the Sellmeier form it asserts the stated error
  against the full parent at its state.
- **Parity:** solprop does not consume these models today, so there is no parity
  fixture.

### 9. Documentation

- `docs/fitting.md`: a section *Compact infrared terms* with the method and the
  tolerance. `docs/data.md` links to it from *Infrared bands (HITRAN)*.
- `docs/gas_models.md`: a table of the compact coefficients per gas, generated by the
  formatter and tested against the stored compact models like the page's other tables.
  It states the qualified intervals, and that the gas absorbs between them.
- A gallery example comparing each compact model with its parent: residuals in sigma
  units over the qualified intervals, and the infrared share of the dispersion.
- Docstrings with units and ranges on `coefficients`, `CompactCoefficients` and the
  loader.

### 10. Compatibility and impact

- Nothing changes for existing users. The defaults, the parent fits, their
  configurations, digests, stored results and robustness records are untouched: the
  reduction has its own configuration (§6).
- `registry.models(gas)` lists one more model per molecule. Code or pages that iterate
  over every model (the gallery's comparisons) will include the compact models. That is
  intended, but it must be reviewed.
- Every refit of a molecular parent requires regenerating its compact model. The
  freshness test fails until that is done, and `AGENTS.md`'s definition of done gains
  the step in Phase 1.
- If OEP 1 is accepted, the compact models carry its identity and capability flags like
  every other model.
- **Sequencing.** OEP 3 changes every molecular parent: its bands, its poles and its
  intervals. Deriving compact models before it would mean deriving them twice. So this
  OEP is implemented after OEP 3's last phase. Nothing here depends on OEP 3's details
  beyond the parents' intervals and representation bound; with one interval and no
  bound, it reduces to the calculation of §8.
- No register entries are superseded.

### 11. Alternatives considered

- **Two poles per band** (Voronin & Zheltikov [1]). It keeps more structure near each
  band, but for CO2 it needs 8 poles for 4 bands and grows with the bands kept. It does
  not reduce across bands, and nothing here needs the band structure.
- **Store more moments per band** (Tang [4], Tang & Poon [5], TC-723). Each band's $M_2$,
  one more column from the same line lists, gives two-sided bounds on the band's
  contribution. An exact two-point rule per band needs four moments of the measure
  $f/\omega^2$: $M_{-2}$, $M_0$, $M_2$ and $M_4$. With $M_2$ alone, a two-point rule
  needs one node fixed, at a band edge (a Radau rule). Either improves the full model
  near the bands, not the compactness. OEP 3 takes this up, as a Gaussian rule of
  several nodes per band built from the cached line lists (§13 Q3).
- **Window polynomials** (Mathar [3]). Valid only inside each window, with no static or
  high-frequency limit and no physical interpretation. Rejected.
- **Line-by-line sums at run time** (Mathar [2]). Exact, but the opposite of compact.
- **Leaving the infrared out of printed formulas** (He et al. [10]). This is the status
  quo that the motivation shows to be wrong by up to 55 % in SF6's GVD.

### 12. Phases

Both phases follow OEP 3's Phase 3, the refitted molecular parents.

1. **The compact models.**
   - `tools/fit/compact.py`, `optmat.gas.models.compact` (the loader) and the registry
     integration.
   - The reduction configurations `data/gas/compact/<gas>.toml` and the stored models
     for CO, CO2, CH4, N2O and SF6.
   - Tests for freshness, tolerance, invariants and reproduction.
   - The `AGENTS.md` regeneration step for the molecular fits.
   - Register entries for the method and tolerance, and for each gas's $K$, bounds and
     achieved residuals.
   - The *Compact infrared terms* section of `docs/fitting.md`.
   - Verify: the gate, and `test-fits` with the re-derivation.
2. **The printed forms.**
   - `coefficients` and `CompactCoefficients` in both forms, and the formatter.
   - The coefficient table on `docs/gas_models.md` (tested) and the gallery example.
   - Verify: the round trip of both forms, the Sellmeier form's stated error against the
     full model at its state, and the gate.

### 13. Open questions

Each question keeps its text as raised. Its disposition follows it (2026-09-30).

1. **The tolerance** (§7): 0.1 sigma in $n - 1$ and $\beta_2$, 0.5 in $\beta_3$. Is that
   the right bar, or should it follow a fraction of the working estimate?

   **Disposition: decided by the author (2026-09-30).** Tighten to 0.1 sigma in
   $\alpha$, $D_2$ and $D_3$ alike. Add 0.01 sigma for rounding, 0.1 of each datum's
   sigma for the reproduction, and exact static values (§7). Every gas meets it with
   margin (worst: CO's $D_3$, 0.057; §8). The working-estimate alternative stays
   rejected, because it is not computable.
2. **CO2's parent.** TC-725 put one of its electronic oscillators on a 202 nm bound;
   TC-726 settled it (the two weak bands below 10 eV as lines). The feasibility numbers
   of §8 were computed against the TC-725 fit and are to be recomputed in Phase 1.

   **Disposition: resolved.** §8 is recomputed against the current parents (after
   TC-726 and TC-730), with their digests. Phase 1 recomputes only if a parent changes
   again, which the freshness test detects.
3. **Per-band $M_2$** (§11): take it up as its own change, before or after this OEP?

   **Disposition: decided by the author (2026-09-30): moved to OEP 3.**
   - A band's higher moments matter at the edges of the infrared windows, which is OEP
     3's subject. There, one pole per band misses the GVD by as much as the whole
     electronic GVD.
   - OEP 3 generalises the question to a Gaussian rule of several nodes per band, built
     from the lines. Its $n = 2$ case is the two-point rule from $M_2$ and $M_4$.
   - The compact models represent whatever the parent stores, so they need no change
     for it.
4. **The printed state** for the Sellmeier form: 273.15 K and 101325 Pa (the
   refractometry convention), or 293.15 K and 100 kPa (the accuracy page's), or both.

   **Disposition: decided by the author (2026-09-30).** Use 293.15 K and 101325 Pa
   (`optmat.conditions.ONE_ATM_20C`). The molecular polarizabilities refer to 293.15 K
   (TC-730), so a table at 273.15 K would put the reference-temperature polarizability
   at another temperature's density without saying so. The export takes any state, so
   this choice fixes only the manual's and the paper's table.
5. **The name.** `compact/<gas>`, or `fit/<gas>/compact`, which keeps the parent in the
   name?

   **Disposition: decided by the author (2026-09-30).** Use `compact/<gas>`, a fourth
   family beside `fit/<gas>` (e.g. `fit/CO2`), `<source>/<formula>` (e.g.
   `old_1971/sellmeier`) and `refractiveindex.info/<gas>/<entry>` (e.g.
   `refractiveindex.info/CO2/Old`).
   - `fit/<gas>` means a stored fit with a configuration and a robustness record.
   - Three-part names are already used by the `refractiveindex.info/<gas>/<entry>`
     comparison models.
   - The parent is named in the record and in `CompactCoefficients`.

The following questions record GPT-6's review of Part I on 2026-09-29. They require
resolution, or an explicit disposition, before acceptance; they do not change the
parent models or select a new physical approximation.

6. **CO2's window contains poles (§1, §3, §7, §8).** The stored parent has 224 fixed
   poles strictly inside 8-12 µm, including a pole at 9.41273 µm with strength
   1.01815e-9. At a relative frequency detuning of 1e-4 that term alone contributes
   about 0.217 au, against an electronic statistical sigma of 1.497e-4 au: about
   1451 sigma. The residual diverges at every parent pole that the compact model
   does not preserve. Thus six poles cannot satisfy a finite error bound at every
   wavelength of that entire interval. A geometric grid, even four times finer,
   cannot establish the stated guarantee. Which pole neighbourhoods are excluded,
   retained exactly, or treated by a separately justified smooth target? Define that
   domain and its provenance before recomputing feasibility. The existing manual's
   broad window claim does not establish a pole-free interval. The snapshot checked
   is `models/fitted/CO2.json`, SHA-256
   `d0be2569a263b73e853cac9e807216736735934cce1acf886a685ed537558fe3`.

   **Disposition: decided by the author (2026-09-30): the windows go to OEP 3, and
   this OEP follows the parents' qualified intervals.**
   - The finding is confirmed on the current parent (§8). The window has no line-free
     interval at all, and the in-window poles contribute more than 0.1 sigma almost
     everywhere in it. The parent itself is not accurate there.
   - OEP 3 defines the parents' windows and their representation:
     - an $n$-node rule per band, from the lines;
     - the weak in-window bands folded into band systems outside;
     - a representation bound measured against the line sum.
     It also takes up the description's and `docs/accuracy.md`'s "valid again between
     the bands, 8-12 µm".
   - The compact models are qualified on whatever intervals the parent qualifies, with
     no pole inside any of them (§4, §7), and weighted by the parent's stated
     uncertainty there. On today's parents that is one interval, 250-1400 nm for CO2.
7. **Representing the domain and constraining new poles (§6, §7).**
   `Model.wavelength_range` and `FittedModel` currently carry one pair of endpoints;
   CO2's stored pair is 250-1400 nm. Copying it warns throughout the proposed second
   interval, while expanding it to 12000 nm declares the intervening bands valid.
   Specify how the compact model records and checks separate qualified intervals
   and Q6's exclusions, and how the export carries them. Bounds spanning all parent
   energies also allow a new effective pole inside a qualified interval. Require a
   pole-free domain check independent of the sampling grid, and constrain energies
   accordingly. Locate the registry integration at
   `optmat.gas.models.registry.models`, the existing entry point; §6 currently names
   an export that `optmat.gas.models` does not provide.

   **Disposition: resolved.**
   - Separate intervals are OEP 3's `Model.intervals`, with `wavelength_range` as their
     envelope. The compact model takes the parent's, and `CompactCoefficients.intervals`
     carries them (§6). Q6's exclusions become the gaps between intervals.
   - Each compact energy is bounded within one absorbing region, between two intervals
     or below the first, so no compact pole can enter a qualified interval. That is
     proved from the bounds, not sampled, and asserted by a test (§4, §8).
   - The registry integration is named in §6 and §7.
8. **Exact constraints in the inner solve (§4, §7, §8).** Plain non-negative least
   squares solves the positivity constraint, but not `sum F = M0` or
   `sum F / Omega**2 = M-2`. Specify an equality-constrained non-negative solve or
   an equivalent parameterisation, its numerical moment tolerance, and the treatment
   of trial energies for which those equalities are infeasible. Penalty rows do not
   make a moment exact. The feasibility table explicitly omits the moment constraints,
   so it does not yet demonstrate the recommended method or its pole counts. Retain
   reproducible inputs, full coefficients, grid and parent digests for the revised
   calculation. When choosing among starts, prefer a candidate meeting every stated
   tolerance: the lowest summed square error can have a worse maximum residual and
   fail a threshold that another candidate passes.

   **Disposition: resolved.**
   - §7 specifies the exact solve (support enumeration with a null-space solve), the
     moment tolerance ($10^{-10}$ relative; achieved $10^{-11}$) and feasibility by the
     bounds.
   - Both moments are now constrained for every gas, which also settles Q10's static
     validators.
   - §8 is recomputed with the constraints and records the parent digests, the grid,
     the starts and the full coefficients.
   - The selection rule prefers candidates that pass.
9. **Absolute sigma and the density convention (§4, §5, §7).** A zero of beta2
   does not imply a zero of its absolute uncertainty. The current SF6 model crosses
   zero at about 1.27122 µm, where the absolute sigma of `D2` is about 0.02835 in
   atomic units, with the parent's usual Birge scaling. Compute sigma directly as
   `sqrt(g.T @ covariance @ g)`; using the relative-uncertainty method at an exact
   zero and then multiplying by the dispersion produces an undefined product.
   Flooring this finite sigma by a fraction of the dispersion elsewhere is not
   justified by the zero crossing and changes the advertised statistical metric.
   State the Birge convention explicitly. Also distinguish the density-independent
   tolerances on alpha and `D2`, `D3` from tolerances on the exact finite-density
   refractivity and dispersion: Lorentz-Lorenz is nonlinear, and the latter need a
   specified state domain and propagated uncertainties, or a dilute-limit label.

   **Disposition: resolved.**
   - §4 now defines the absolute sigma from the covariance with the Birge convention
     stated, and drops the floor. SF6's value is confirmed: 1271.2 nm, 0.0284 au.
   - The tolerance is labelled dilute-limit.
   - The finite-density effect is bounded analytically to first order and measured
     ($8\times10^{-4}$ of itself at 10 bar), and a test repeats it at 10 bar.
10. **Observable-specific checks and static validators (§5, §7, §8).** A pointwise
    error in units of model sigma does not define a tolerance on a dataset's
    chi-squared or on a phase-matching ratio. For fixed data-normalised residuals
    `r` and prediction changes `delta`, the chi-squared change is
    `2 * r.T @ delta + delta.T @ delta`; it depends on the original residuals and
    the dataset size. Define the allowed changes per observable, the handling of
    dataset scales and partner models, and tests at the actual harmonic wavelengths.
    The parent's validators also include static dielectric values for CH4 and a
    total static polarizability for CO. These depend on `M-2` even though those
    compact models are qualified only at optical wavelengths. Either constrain or
    separately test that moment, or explicitly exclude those extrapolation checks
    from the preservation promise. Marchetti's CO2 data and several Shelton ratios
    already enter the parent fits: they are useful reproduction checks, but are not
    independent validation of the parent's physical accuracy.

    **Disposition: resolved.**
    - The reproduction tolerance is per datum: 0.1 of its own sigma, inside the declared
      range. The fitted scales and partner models stay fixed, and the checks run at the
      data's own wavelengths, harmonics included. Chi-squared changes are reported, not
      bounded (§7).
    - $M_{-2}$ is constrained, so the static validators are exact (§8).
    - §5 now calls Marchetti's values and Shelton's ratios reproduction checks.
11. **Configuration digest and regeneration (§6, §10, §12).** `load_config` hashes
    the complete TOML bytes. Adding `[compact]`, or later changing its grid or pole
    count, changes the digest in both the stored parent and its robustness record;
    `tests/test_stored_fits.py` checks both against the current configuration.
    Therefore §10's promise to leave the parent records untouched is incompatible
    with the proposed storage without a further design choice. Prefer separate
    reduction configuration, or specify an explicit digest migration and which
    records must be regenerated. Put the regeneration rule in Phase 1, when the
    compact records and freshness checks first become part of the repository.

    **Disposition: resolved.** The reduction configuration is a separate file,
    `data/gas/compact/<gas>.toml`, whose digest the compact record carries. The parent's
    configuration, digest and records are untouched (§6, §10). The regeneration rule and
    the freshness test are in Phase 1 (§12).
12. **The quadrature construction (§4, §7, §11).** State the quadrature variable and
    measure. Gaussian quadrature in `x = omega_b**2` with weights `f_b` matches the
    even moments governing the optical expansion; quadrature in `omega_b` matches
    a different sequence. Grouping bands and preserving each group's `M0` and
    `M-2` is not that Gaussian construction. An independent six-node Gaussian
    calculation in squared frequency on the current CO2 parent does fail the optical
    tolerance (on 1600 geometric points, about 24 sigma in `D2` and 7266 in `D3`),
    supporting the motivation for an optical fit, but not the claimed equivalence
    of the two trials. In §11, the three moments `M-2`, `M0`, `M2` do not determine
    a unique unrestricted two-pole rule: identify the fourth moment, or a fixed
    endpoint and the corresponding quadrature construction, before calling it an
    exact two-point quadrature.

    **Disposition: resolved.**
    - §4 states the variable ($x = \omega_b^2$) and the measure ($f_b$), and the moments
      matched. It drops the grouped-band trial, and quotes the Gaussian rule's failure,
      reproduced independently: 23 sigma in $D_2$ and 7320 in $D_3$ on the current
      parent.
    - §11 now says that an exact two-point rule needs $M_4$ as well, or one fixed node
      (Radau).
13. **Qualifying the printed coefficients (§6, §8, §9).** Specify the return type,
    coefficient ordering and units, qualified intervals, number density used at the
    printed state, parent identity and error metadata. A printed table needs enough
    significant digits to preserve the error bound after its decimal strings are
    read back; testing the stored binary coefficients alone does not test the main
    publication use case. Round-trip both printed forms and budget rounding error
    together with reduction error. Define whether the Sellmeier form's stated error
    is an absolute or relative error, over which wavelengths and for which
    observables; test the complete difference against the parent, including the
    omitted Lorentz-Lorenz and virial contributions. The leading `(n-1)/6` estimate
    is an expansion, not a bound over every state allowed by an equation of state.

    **Disposition: resolved.**
    - §6 specifies `CompactCoefficients`: its fields, ordering, units, qualified range,
      number density, parent identity and error metadata.
    - The rounding budget is 0.01 sigma, and 10 significant digits meet it (CO2 needs
      them; §8).
    - Both forms are round-tripped from their decimal strings.
    - The Sellmeier form's error is the largest relative difference in $n - 1$ and
      $\beta_2$ from the full parent at the stated state over the qualified intervals,
      computed directly rather than from the $(n-1)/6$ estimate.
14. **Reconcile the feasibility claims with the parent (§1, §3, §8, Q2).** Q2
    acknowledges the old CO2 fit, but the abstract still presents its pole counts
    and errors as the proposed method's result. The current CO2 model has four
    leading electronic terms, so six infrared terms make ten terms, not the stated
    maximum of eight. Its fixed poles reach about 6996.69 cm^-1, beyond the
    6972 cm^-1 used in the edge argument. Update the parent inventory, separate
    provisional unconstrained results from qualified results, and recompute after
    Q6 and Q8 are settled. The present review could not reproduce the original
    feasibility table because its full coefficients and calculation are not recorded
    in the proposal or a public derivation tool.

    **Disposition: resolved.**
    - §3 has the current inventory: electronic terms, fixed poles, the highest pole and
      the range.
    - The abstract and §8 give the constrained results on the current parents, with 4
      to 10 terms in all. The edge argument uses 6996.7 cm⁻¹.
    - The unconstrained trial table is withdrawn.
    - §8 records the full coefficients and digests.

### 14. References

1. A. A. Voronin and A. M. Zheltikov, "The generalized Sellmeier equation for air,"
   *Sci. Rep.* **7**, 46111 (2017). doi:10.1038/srep46111 (`voronin_2017`)
2. R. J. Mathar, "Calculated refractivity of water vapor and moist air in the atmospheric
   window at 10 µm," *Appl. Opt.* **43**, 928 (2004). doi:10.1364/AO.43.000928
   (`mathar_2004`)
3. R. J. Mathar, "Refractive index of humid air in the infrared: model fits," *J. Opt. A:
   Pure Appl. Opt.* **9**, 470 (2007). doi:10.1088/1464-4258/9/5/008 (`mathar_2007`)
4. K. T. Tang, "One-term Sellmeier formula for dispersion of dilute gases," *J. Opt. Soc.
   Am.* **62**, 644 (1972). doi:10.1364/JOSA.62.000644 (`tang_1972`)
5. K. T. Tang and K. K. Poon, *J. Opt. Soc. Am.* **64**, 1582 (1974).
   doi:10.1364/JOSA.64.001582 (`tang_1974`)
6. I. E. Gordon, L. S. Rothman, R. J. Hargreaves et al., "The HITRAN2024 molecular
   spectroscopic database," *J. Quant. Spectrosc. Radiat. Transfer* **353**, 109807
   (2026). doi:10.1016/j.jqsrt.2026.109807. A data service, not a curated source record;
   the provenance of the band tables is `data/external/hitran/PROVENANCE.md`.
7. S. Marchetti and R. Simili, *Infrared Phys. Technol.* **50**, 78 (2007).
   doi:10.1016/j.infrared.2006.07.003 (`marchetti_2007`)
8. D. P. Shelton, *Phys. Rev. A* **42**, 2578 (1990). doi:10.1103/PhysRevA.42.2578
   (`shelton_1990`)
9. D. P. Shelton and V. Mizrahi, *Chem. Phys. Lett.* **120**, 318 (1985).
   doi:10.1016/0009-2614(85)87065-2 (`shelton_1985_sf6`); D. P. Shelton, *J. Chem.
   Phys.* **85**, 4234 (1986). doi:10.1063/1.451794 (`shelton_1986_co2`)
10. Q. He, Z. Fang, O. Shoshanim, S. S. Brown and Y. Rudich, *Atmos. Chem. Phys.* **21**,
    14927 (2021). doi:10.5194/acp-21-14927-2021 (`he_2021`)

## Part II — Implementation record

### Before Phase 1 — the scope on OEP 3's windowed parents

- **Finding (TC-798):** OEP 3's Phase 3 gives the molecular parents windows beyond their
  optical ranges. There the fixed infrared terms carry no statistical uncertainty, so
  §7's tolerance, a tenth of the parent's sigma, would demand the band rules themselves.
  A feasibility calculation on the refitted CO, with §7's method:
  - five poles per absorbing region miss by 56 sigma in $D_2$ and 2300 in $D_3$;
  - with the representation bound as the tolerance in the windows, about 2 sigma at
    best;
  - on the optical range alone, three poles meet the tolerance in $\alpha$ and $D_2$.
- **Decision (the author, 2026-10-01):** a deviation from §4 and §6, recorded here as
  OEP 0 requires.
  - The compact models are qualified on the parents' optical intervals only: their
    qualified intervals inside the optical range (`range_nm`), with §7's tolerance.
  - The windows beyond are left to the full model, and a compact model's range warning
    names the parent for them.
  - Compact poles are kept outside every optical interval; they may lie in the parents'
    windows.

### Phase 1 — The compact models

- **Scope:** Phase 1 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0002-compact`.
- **Delivered:**
  - {mod}`optmat.gas.models.compact`: the reduction (`derive`), the residuals, the
    reproduction check and the loader, `compact_models()`, joined to the registry;
  - `tools/fit/compact.py`;
  - the configurations `data/gas/compact/<gas>.toml` and the stored models for CO,
    CO2, CH4, N2O and SF6;
  - `FittedModel.dispersion_sigma_au`, the absolute sigma of §4;
  - tests of freshness, tolerance, invariants, reproduction and the finite density, and
    the re-derivation in the `fits` lane;
  - the section *Compact infrared terms* of {doc}`../../fitting`, with a tested table;
  - the regeneration step in `AGENTS.md`.
- **Deviations from the proposal:**
  - The intervals are the parents' optical ones (Part II, *Before Phase 1*). CH4's and
    N2O's are shortened further, to 240-1070 nm and 320-1500 nm, as decided with the
    author (below; TC-800).
  - The absorbing regions include one above the intervals: OEP 3's fold puts a share of
    each top band system on a node beyond the parent's short-wavelength end, at 107-148
    nm. No compact model gives it a pole; the other poles take its share within both
    moments.
  - A region given as many poles as it holds keeps the parent's poles unchanged; SF6's
    three are carried through this way, reordered by energy.
  - The absolute sigma of §4 is a new public method, `FittedModel.dispersion_sigma_au`;
    `dispersion_sigma_relative` now divides it by the dispersion.
- **Decisions made:**
  - Poles per region: CO 3, CO2 5, CH4 5, N2O 5 below the intervals; none above.
  - Starting points: §7's, with four seeded random starts (twelve for CO2).
  - The compact model's status is `comparison`, with the reason that it is derived for
    print and never a default.
  - The reproduction counts a static value apart from a dataset's dynamic ones.
- **Results (TC-799):**

  | Model | Interval | Compact poles | $\alpha$ | $D_2$ | $D_3$ |
  |---|---|---:|---:|---:|---:|
  | `compact/CO` | 167-2000 nm | 3 | 0.0017 | 0.0057 | 0.049 |
  | `compact/CO2` | 250-1400 nm | 5 | 0.0022 | 0.0015 | 0.046 |
  | `compact/CH4` | 240-1070 nm | 5 | 0.0007 | 0.0006 | 0.010 |
  | `compact/N2O` | 320-1500 nm | 5 | 0.0006 | 0.0002 | 0.069 |
  | `compact/SF6` | 230-5000 nm | 3 (kept) | 0 | 0 | 0 |

  - One pole fewer fails in every reduced case.
  - $M_0$ and $M_{-2}$ hold to $10^{-15}$.
  - Every datum inside the intervals moves by at most 0.003 of its sigma, and the static
    values are exact to $6\times10^{-16}$.
  - At 10 bar the finite-density correction is $10^{-3}$ of the reduction error.
- **Register entries:** TC-799 (the method and results), TC-800 (the finding on CH4 and
  N2O).
- **Verification:**
  - The full gate, and the `fits` lane's re-derivation of all five.
  - The golden values of every existing model are unchanged; the compact models are added.
- **Sources checked:** none new; the method's references are curated (§14).
- **Limitations:**
  - The compact models do not cover the parents' infrared windows (Part II), nor, for
    CH4 and N2O, the near infrared beyond 1070 and 1500 nm (Q15).
  - The printed forms, `coefficients` and the coefficient table are Phase 2.
- **Decided with the author: CH4's and N2O's compact intervals (TC-800).** Their
  optical ranges are split by narrow absorbing gaps, whose band-rule poles lie within
  0.2 to 3 cm⁻¹ of the interval edges; there a single parent pole contributes up to
  $10^{16}$ of the parent's sigma in $D_3$, and no few poles follow it. The options
  were:
  - end the compact intervals at 1070 nm (CH4, 5 poles, worst 0.010 sigma) and 1500 nm
    (N2O, 5 poles, 0.069), leaving the near infrared beyond to the full models;
  - keep the parents' gap poles unchanged, at 36 (CH4) and 20 (N2O) poles and more
    beyond the last interval: hundreds of terms, not a printable model;
  - leave CH4 and N2O without a compact model.

  *Decided (the author, 2026-10-01):* the shortened intervals.

### Phase 2 — The printed forms

- **Scope:** Phase 2 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0002-compact`.
- **Delivered:**
  - `coefficients`, `CompactCoefficients` and `format_coefficients` in
    {mod}`optmat.gas.models.compact`, in both forms;
  - the section *The compact coefficients* of {doc}`../../gas_models`, both forms for each
    molecule, tested against the formatter;
  - the gallery example *The compact models against their parents*.
- **Deviations from the proposal:**
  - Eleven significant digits, not ten: on the present parents ten change CO2's $D_3$
    by 0.0101 of its sigma, above the 0.01 budget; eleven give at most 0.0017.
  - The Sellmeier form's $\beta_2$ error is stated relative to the largest $|\beta_2|$
    over the intervals: $\beta_2$ changes sign inside the intervals of CO, CO2, N2O and
    SF6, where a pointwise relative difference is unbounded.
  - The terms' kinds include `rotational`, OEP 3's rotational term.
- **Decisions made:**
  - The printed covariance's parameters are the electronic terms' $f_k$ (or $B_k$) in
    printed order, then the effective terms' $\lambda_k$ (or $C_k$), Birge-scaled as the
    parent's.
  - The page prints the Sellmeier form at 101325 Pa and 293.15 K (Q4).
- **Results (TC-801):** the Sellmeier form's error against the full parent at that state
  is 8.3e-5 to 2.1e-4 in $n - 1$ and 1.7e-4 to 2.6e-4 of the largest $\beta_2$; the
  printed covariance reproduces the parent's sigma to rounding.
- **Register entries:** TC-801.
- **Verification:** the full gate; the round trip of both forms from their decimal
  strings; the Sellmeier form's stated error recomputed from the printed values.
- **Sources checked:** none new.
- **Limitations:** as Phase 1's. With both phases implemented, `Final` awaits the
  author.

## Part III — History

- 2026-09-29: Drafted.
- 2026-09-29: Q2 updated: CO2's parent settled by TC-726.
- 2026-09-29: GPT-6 (AI) reviewed all of Part I against the current models,
  configurations, fitting and uncertainty code, tests, manual and register. Recorded
  unresolved findings as Q6-Q14; the recommendation is to revise before acceptance.
  Numerical checks counted CO2's poles inside the proposed window, evaluated their
  contribution relative to the covariance, propagated SF6's uncertainty at its
  zero-dispersion wavelength, and tested Gaussian quadrature in squared frequency.
  Primary-source pages checked: Tang 1972, p. 645 (the one-term bound); Mathar 2007,
  p. 473 (the instability of state-dependent Sellmeier parameters); Voronin and
  Zheltikov 2017, pp. 2-3 (the dilute oscillator model and transparency restriction).
  HITRAN's official citation page confirms reference [6]; its DOI is
  [10.1016/j.jqsrt.2026.109807](https://doi.org/10.1016/j.jqsrt.2026.109807).
  SciPy's [nnls documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html)
  confirms that its constraints are non-negativity only. Peck's original, the
  measurement papers, and the HITRAN line lists were not freshly checked; no
  independent reproduction of the original feasibility table is claimed. No data,
  fit or physical convention was changed.
- 2026-09-30: Revised by the drafting agent (Opus 5.5) in response to GPT-6's review.
  Every question now has a disposition. Q2 and Q7-Q14 are resolved in the text. Q1, Q3,
  Q4 and Q5, and the scope narrowing of Q6, are recommendations for the author to
  confirm. The changes:
  - Scope: the compact models are qualified on each parent's declared range only, and
    CO2's 8-12 µm window is out.
  - Method: both $M_0$ and $M_{-2}$ are exact for every gas, through an exactly
    constrained inner solve. Sigmas are absolute from the covariance, Birge-scaled, with
    no floor.
  - Tolerance: 0.1 sigma in $\alpha$, $D_2$ and $D_3$, plus rounding, per-datum and
    static tolerances.
  - Storage and interface: a separate reduction configuration, and a specified
    `CompactCoefficients`.
  - Recomputed against the parents stored after TC-730, with digests and full
    coefficients recorded. The unconstrained 2026-09-29 table is withdrawn.
  - Also checked on the current parents: CO2's window, per datum for every parent
    dataset and validator, rounding, finite density and the Sellmeier form.
  - Peck's reference, which was not curated and served only as the format Mathar
    compares with, is dropped.
  - The sigmas were checked against the package's own methods. GPT-6's SF6 and
    Gaussian-quadrature numbers were reproduced.
  - This revision is not an independent review: the author's human review is still
    required for `Reviewed`. No data, fit or physical convention was changed.
- 2026-09-30: The author decided Q1, Q3, Q4, Q5 and Q6.
  - Q3 and CO2's window move to the new OEP 3 (the molecular gases between their
    infrared bands), which is to be implemented first.
  - This OEP now requires OEP 3, and is written to follow the parents' qualified
    intervals: compact poles are bounded to the absorbing regions between them, the
    weights include OEP 3's representation uncertainty, and `CompactCoefficients`
    carries the intervals.
  - The feasibility numbers are labelled as those of today's single-interval parents,
    to be recomputed in Phase 1.
  - The author's decisions on these questions are not a full review of Part I, and the
    review status is unchanged.
- 2026-09-30: Aligned with OEP 3's revision for its review. OEP 3's representation
  uncertainty is now a representation bound, which is stated beside the statistical
  sigma and not combined with it in quadrature (§4, §7). The weights are the parents'
  statistical sigma. OEP 3's rotational band system is carried outside the reduction. No
  number changes on today's parents, where the bound is absent.
- 2026-09-30: John Travers (human) reviewed Part I. With GPT-6's review, the review
  status is `Reviewed`.
- 2026-09-30: Accepted by John Travers. Part I is frozen.
- 2026-10-01: The author decided the scope on OEP 3's windowed parents before Phase 1:
  the optical intervals only (Part II; TC-798).
- 2026-10-01: Phase 1 implemented (TC-799). The author decided CH4's and N2O's compact
  intervals: 240-1070 nm and 320-1500 nm (Part II; TC-800).
- 2026-10-01: Phase 2 implemented (TC-801).
- 2026-10-01: The author checked the compact models and their tables. Final.
