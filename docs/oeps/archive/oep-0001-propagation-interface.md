# OEP 1 — The propagation-code interface: identity, overrides and the data for solprop's SEP 3

| | |
|---|---|
| **OEP** | 1 |
| **Title** | The propagation-code interface: identity, overrides and the data for solprop's SEP 3 |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | GPT-6 (AI); Opus 5.5 (AI); John Travers (human) |
| **Type** | Feature |
| **Created** | 2026-09-29 |
| **Resolution** | 2026-09-30: accepted by the author, John Travers, after the human and AI reviews; merged directly to `main` |
| **Requires** | — |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phases 1-5 merged (`feature/oep-0001-parity`, `-propagation-page`, `-release`, `feature/oep-0001-identity`, `-overrides`, `-derivative`, `-eos`, `-transport`, `-air`, `-raman-groups`, `-rate-energies`, `-silica-bundle`, `-solid-tables`) |
| **Manual pages** | {doc}`../../propagation`; {doc}`../../identity`; {doc}`../../transport`; {doc}`../../models` (the air models) |

## Part I — Proposal

### 1. Abstract

solprop's enhancement proposal SEP 3 makes optmat the only source of the propagation
code's built-in material data:

- the gases' dispersion, their Kerr response, and the equation of state and transport;
- the solids' index;
- the gas and silica Raman responses;
- the gas ionisation rates.

For that to be possible, optmat needs a stable, documented surface that a propagation
code can build on. This OEP proposes it:

- a declared public surface with a snapshot test;
- an identity on every model that changes whenever the model's numbers can;
- machine-readable capability flags;
- construction from explicit parameters, and one-line overrides;
- the SI derivative of every gas polarizability;
- an inverse of the equation of state on a declared phase domain, and sourced
  transport properties;
- named air compositions and an air model;
- a group view of the Raman responses;
- the energy lost per ionisation event of each rate, channel by channel;
- a consistent silica nonlinear bundle;
- a refusal option on the solid tables;
- re-recorded parity fixtures, a first tagged release and a manual page for propagation
  codes.

Success is judged in three ways. solprop's SEP 3 adapters can be written against the
declared surface alone. Every new quantity is verified against an analytic limit or an
independent reference. Nothing optmat already validated changes its numbers except where
§10 says so.

The review of 2026-09-30 checked every open question against the code at `ea0dcaf` and
ran the checks that could be run then: the paired silica gain (§4.7), the equation of
state across phases (§4.2), the low-density widths (§4.5) and the rate channels (§4.6).
Its resolutions are written into the sections below, and §13 records each question's
outcome.

### 2. Motivation and scope

solprop's SEP 3 ([solprop SEP 3][sep3]) deletes solprop's own material tables and
ionisation formulas and takes everything from optmat. Its rules that bear on optmat are:

- **D1:** solprop keeps no rate formulas. Every rate, for built-in or custom species,
  comes from an optmat rate model.
- **D2:** Raman responses enter solprop's existing kernel as one term per damping group
  and tensor family.
- **D4:** silica's scalar Raman default is Hollenbeck and Cantrell's shape, with the
  Raman fraction of optmat's Kerr fit, $f_R = 0.219$ (TC-696). The fraction must
  therefore be available on its own.
- **D5:** solprop keeps linear density scaling, $\chi = N\alpha/\varepsilon_0$, with
  optmat's $\alpha$ and a declared, measured error.
- **D6:** larger optmat design decisions are reviewed and implemented in optmat
  beforehand. This OEP is that design.
- **D7:** quick experiments stay one line. Every model solprop uses accepts explicit
  parameter overrides, and the result still knows what it is.
- **D8 to D11:** the energy-loss convention, the air model, the weak rate defaults and
  the silica tensor pairing, for which §7 proposes optmat's side.

The use case is concrete. A propagation run must be able to:

- say exactly which material models it used, and refuse to resume with different ones;
- warn when a model is provisional;
- let a user try "this rate with a different structure coefficient" or "these Raman
  widths times 0.8" without editing optmat.

**In scope** is everything in the abstract, keyed to solprop SEP 3's upstream gaps
G1–G10 (its §3.2) in §3.

**Out of scope** is anything solprop can do itself:

- tabulation, and rate tables;
- the kernel mapping of damping groups;
- mean free path and sound-speed formulas;
- the vacuum limit;
- solprop's $n_0$ conventions;
- the reference-temperature rule.

Also out of scope:

- **The mixture law of the gas index at finite density** (solprop's SEP 11). The virial
  mixture law with its unlike-pair coefficients, and the density of a mixture, are
  [OEP 4](oep-0004-mixture-index.md), which builds on §6.4's derivative and §6.5's
  equation-of-state object.
- **New gases and new Raman physics.** Line mixing, N2O's rotations and the other open
  gaps listed on the {doc}`roadmap <../../roadmap>` are not addressed here.
- **The molecules' infrared windows**, and the representation of their bands
  ([OEP 3](oep-0003-infrared-windows.md)). This OEP declares a model's validity as a set
  of qualified intervals on the surface (§4.1), with one interval per model. OEP 3
  gives the molecular fits several.

Two constraints hold throughout:

- **No dependence on solprop.** Nothing here makes optmat depend on solprop, or adds a
  code path that only solprop uses.
- **Conventions unchanged.** optmat's quantities stay SI and per particle. Helpers
  already written in solprop's convention, such as `gamma3`, stay documented as such.

### 3. Context

optmat at commit `d57963b`, re-checked at `ea0dcaf`, against solprop SEP 3's gaps:

| Gap | optmat today | Section |
|---|---|---|
| Public surface | solprop's adapters will import from about a dozen modules (listed below). Nothing marks which names are contract, and several useful builders are private | §6.1 |
| G4, G9. Identity and flags | Gas fits carry a configuration SHA-256 in their citation string. It hashes the fit's TOML bytes only, not the datasets, absorption spectra or band tables the fit used, nor the stored terms that evaluation reads. Solid fits carry `config_sha256` in their JSON without exposing it. `KerrModel`, `SilicaRamanResponse`, `GasRamanResponse`, `IonisationRateModel` and `MolecularRateModel` have no identity. Provisional and untested statuses (N2O's Raman bands, SF6's rate, the suggested α of Kr and Xe) exist only in prose. The models are frozen dataclasses, but several hold writeable arrays or dicts (`FittedModel.covariance`, `FormulaModel.terms`, `MeasuredAbsorption`'s grids, every response's arrays, the rate models' `sources`), and the registries return shared cached instances | §6.2 |
| G5. Explicit parameters | The rate models are frozen dataclasses with public fields, but they are reachable only through `ionisation_rate_model(gas, λ, model)` and `molecular_rate_model(...)`, typed on the `Gas` enum, so a custom species cannot be expressed. `IonisationRateModel` stores no derived value; `MolecularRateModel` stores two (`structure_factors`, and `coefficient_squared` for the PPT and ADK forms), which `dataclasses.replace` would leave stale. `GasRamanResponse` is built only from curated data. The solids' `FormulaModel` and `TabulatedModel` exist, but not as documented constructors | §6.3 |
| G1. Polarizability derivative | `complex_polarizability(λ)` exists in SI (C m²/V) on every gas model. Only `OscillatorModel.dpolarizability_domega_au` gives a derivative, real and in atomic units. It is analytic except for O3's measured absorption, whose `real_derivative_au` is a repeated `numpy.gradient` of the grid, linearly interpolated, and not the derivative of the evaluated spline (§4.1) | §6.4 |
| G2. EOS inverse and transport | `optmat.gas.eos.number_density(gas, State)` and `compressibility` exist, through pyfluids (CoolProp 8.0.0). `State` needs a positive pressure. There is no inverse and no transport property. A (p, T) input past saturation silently returns the liquid's density: SF6 at 25 bar and 294 K gives Z = 0.107, CoolProp phase `Liquid` | §6.5 |
| G3. Compositions | `optmat.gas.composition.composition("dry air")` is CIPM-2007 dry air (Picard et al. [1]): N2 0.780848, O2 0.20939, Ar 0.009332, CO2 0.0004 and eight traces, summing to 1 − 2.9 × 10⁻⁶. `raman_response("dry air")` keeps N2, O2 and Ar and renormalises them. solprop's two recipes are NOAA's rounded fractions, renormalised: N2/O2/Ar, and N2/O2, which its interface string `"air"` uses | §6.6 |
| G6. Raman groups | `GasRamanResponse` carries its lines' species, band, branch, J, frequency, isotropic and anisotropic amplitudes, width law, Doppler width, shift, damping group, and density and temperature ranges ({doc}`../../raman_response`, TC-685 to TC-692). It has no line identifier. A density ≤ 0 raises, and a density below the range warns | §6.7 |
| G7. Rate energies | `IonisationRateModel.channels` holds each channel's ionisation energy and weight. `MolecularRateModel` has `ionization_energy`, its kind (vertical in the default) and `adiabatic_ionization_energy`. `maximum_field`, `uncertainty_factor` and `keldysh` exist ({doc}`../../ionisation_rates`, {doc}`../../molecular_ionisation_rates`, TC-652, TC-677 to TC-680). No per-channel rate is public, and the default "consistent" averaging factor is applied to the channel sum. The rates are stated for linear polarisation only | §6.8 |
| G8. Silica bundle | `kerr/SiO2` gives $n_2(\lambda)$, $n_{2,e}(\lambda)$ and $f_R = 0.219 \pm 0.026$ (TC-696). `silica_raman_response("hollenbeck_2002" \| "lin_2006" \| "blow_1989")` exists, with Hollenbeck paired to the Kerr fit. `raman_coefficient` and `electronic_susceptibility` include $n_0^2$ (TC-692) ({doc}`../../silica`) | §6.9 |
| G4, G9 for solids | Every solid solprop ships has a recommended model with exact derivatives: `fit/SiO2`, `daimon_2002/sellmeier_20C`, `dodge_1984` o and e, `malitson_1964`, `weber_2003` o and e, `schott_2023`, and the tables `table/Ag` and `table/Al`. `default_model` refuses a missing axis. The metal tables hold their end values with a warning, and N-BK7 is lossless ({doc}`../../crystals`, {doc}`../../metals`) | §6.10 |
| Parity fixtures | `tests/data/solprop_*.json` record solprop's current values: gamma3, sellmeier, raman, raman_response, ionisation, molecular_rates and solids (TC-439, TC-690) | §6.11 |
| G10. Distribution | The version comes from hatch-vcs. There are no tags and no PyPI release, and the changelog has only *Unreleased* | §6.12 |

The modules solprop's adapters will import are:

- `optmat.gas.models.registry`, `optmat.gas.nonlinear` and `optmat.gas.eos`;
- `optmat.gas.ionisation`, `optmat.gas.ionisation_rates` and
  `optmat.gas.molecular_ionisation_rates`;
- `optmat.gas.raman_response` and `optmat.gas.composition`;
- `optmat.solid.models.registry`, `optmat.solid.kerr` and `optmat.solid.raman_response`;
- `optmat.exceptions` and `optmat.conditions`.

The {doc}`roadmap <../../roadmap>` item *Raman and ionisation parameters for propagation*
(TC-439) and its *solprop* section describe the integration this OEP prepares. The
species typing of datasets is to be generalised under TC-317.

### 4. Theory

#### 4.1 The polarizability derivative

A gas model gives the complex dynamic polarizability $\alpha(\omega)$ per particle, in
C m²/V, at vacuum wavelength $\lambda = 2\pi c/\omega$. The derivative wanted is
$\partial\alpha/\partial\omega$ at fixed state, in C m² V⁻¹ s. It enters the group index
and every frame derivative of a propagation code. With $\omega = 2\pi c/\lambda$,
$d\lambda/d\omega = -\lambda/\omega$, so

$$
\frac{\partial\alpha}{\partial\omega} = -\frac{\lambda}{\omega}\,\frac{\partial\alpha}{\partial\lambda}.
$$

**Scope.** The derivative is provided on every model the gas registry returns
(`default_model`, and every named fitted, published and database model: 106 at
`ea0dcaf`). `optmat.gas.euv.ComplexModel` is not a registry model: it has no range and
no `complex_polarizability`, and it stays out of scope (§13, Q13).

The model classes are as follows.

- **The oscillator and fitted models** ({doc}`../../models`): in atomic units,
  $\alpha(\omega) = \sum_i f_i/(E_i^2 - \omega^2) + \alpha_\text{off}
  + \alpha_\text{ret}\,\omega^2 + \alpha_\text{meas}(\omega)$, where $\alpha_\text{ret}$
  is helium's retardation term (TC-728) and zero elsewhere, and the fixed infrared bands
  of the molecular fits are terms of the sum, however many poles each band has (one
  today, a few from OEP 3). So
  $\partial_\omega\alpha = \sum_i 2 f_i\,\omega/(E_i^2-\omega^2)^2
  + 2\alpha_\text{ret}\,\omega + \partial_\omega\alpha_\text{meas}$.
  - The conversion to SI multiplies $\alpha$ by the atomic unit of polarizability,
    $e^2 a_0^2/E_h$.
  - The conversion multiplies $\partial_\omega$ by $\hbar/E_h$ (from the hartree energy
    to rad/s).
- **The measured absorption** (O3's only, {doc}`../../gas_models`). The derivative is that
  of the evaluator actually used, branch by branch:
  - inside the grid, the real part is a cubic spline of the Kramers–Kronig transform, and
    its derivative is the spline's own;
  - outside the grid, the real part is a trapezoid sum over the grid, and its derivative
    is the same sum of the differentiated integrand, $2\nu\,\sigma(\nu')/(\nu'^2-\nu^2)^2$;
  - the imaginary part is $\varepsilon_0 c\,\sigma(\nu)/\omega$ with $\sigma$ linearly
    interpolated, so its derivative uses the slope of the segment. At a knot, where the
    interpolant has a kink, the convention is the mean of the two slopes; at the grid's
    upper end the imaginary part drops to zero, and the derivative there is that of the
    segment below.

  The existing `real_derivative_au`, a repeated `numpy.gradient` of the grid, differs
  from the spline's derivative by up to 4.6 × 10⁻⁵ of its maximum inside O3's range, has
  the wrong sign outside the grid, and is not zero at $\omega = 0$ (TC-732). It is not
  reused. The
  grid ends at 96 787 cm⁻¹ (103 nm), outside O3's range (185–3800 nm), where the real
  part is discontinuous; there the derivative is that of the branch evaluated, and the
  range warning applies.
- **The published formulas and database entries** ({doc}`../../gas_models`): each formula is
  a closed form in $\lambda$ or in wavenumber. Its refractivity is reduced to a
  polarizability by inverse Lorentz–Lorenz with the refractivity virial, at the number
  density of the formula's reference state. That density is a constant of the model, so
  the derivative is the formula's analytic derivative times the Jacobian of the
  reduction. With $\alpha = 3\varepsilon_0\chi/[N P(\chi+3)]$, where $P$ is the pair
  factor of the refractivity virial at that density (`optmat.gas.lorentz_lorenz`),
  $\partial\alpha/\partial\chi = 9\varepsilon_0/[N P(\chi+3)^2]$ and
  $\partial\chi/\partial(n-1) = 2n$.

No registry model needs a numerical derivative. The complex step
$f'(x) \approx \operatorname{Im} f(x+ih)/h$ would not apply in any case: $\alpha$ is
complex-valued.

**Range semantics** are those of `complex_polarizability`: the same `ModelRangeWarning`
outside every qualified interval (below). The derivative is computed in $\omega$, not through the chain
rule, which is $0\cdot\infty$ at $\lambda = \infty$. In the static limit it is exactly 0
for every class, since the real parts are even in $\omega$ and the imaginary part
vanishes there; O3's spline, whose derivative at $\omega = 0$ is 1.6 × 10⁻⁹ a.u. from
its end condition, is given 0 there by definition. The
seven database entries of form 2 that return NaN at $\lambda = \infty$ (Borzsonyi's and
Ermolov's refractiveindex.info formulas) keep that behaviour.

**Qualified intervals.** A gas dispersion model's validity is a tuple of disjoint,
ascending intervals of vacuum wavelength, `Model.intervals`, and a wavelength outside
every one of them warns.
- Every model has one interval in this OEP, its present `wavelength_range`, so no
  warning changes.
- `wavelength_range` is kept as the envelope, from the first interval's lower end to
  the last one's upper end, and documented as such.
- OEP 3 gives the molecular fits several intervals, with gaps between them where a band
  absorbs or its representation is not accurate. There the envelope is not the
  validity, and a wavelength between two intervals warns with the gap's reason.

Declaring the tuple now keeps that later change within the surface, rather than a change
of it. A consumer tests validity with `intervals` or through the warning, never against
`wavelength_range` alone.

#### 4.2 The equation of state and its inverse

optmat's number density is $N = p/(Z(p,T)\,k_B T)$, with the compressibility factor $Z$
from the real-gas equation of state that `optmat.gas.eos` uses. The inverse wanted is
the pressure at given $N$ and $T$. It solves

$$
F(p) = p - N k_B T\, Z(p, T) = 0.
$$

**The phase domain.** With a (p, T) input, CoolProp returns the stable phase, and the
forward law today silently returns a liquid's density: SF6 at 294 K gives
$Z = 0.688$ (`Gas`) at 20 bar and $Z = 0.107$ (`Liquid`) at 25 bar, and CO2 at 280 K and
50 bar is liquid. optmat is a library of gases, so the forward law and its inverse share
one declared domain: the gas below saturation ($T < T_c$, $p < p_\sigma(T)$) and the
supercritical fluid ($T \ge T_c$). A state CoolProp reports as `Liquid`,
`SupercriticalLiquid` ($T < T_c$, $p > p_c$) or two-phase raises `StateError`, forward and
inverse alike. This changes the forward law only where it now returns a condensed phase
(§10, TC-731).

**The solver.** On that domain the density increases strictly with pressure at fixed
temperature, because $(\partial p/\partial\rho)_T > 0$ in any mechanically stable single
phase. $N(p)$ is therefore monotonic, and the root of $F$ is unique on a bracket:

- below $T_c$, the bracket is $(0, p_\sigma(T))$, and a density at or above the saturated
  vapour's, $N \ge N(p_\sigma, T)$, is refused as condensed;
- at or above $T_c$, the bracket is $(0, p_\text{max})$, the equation of state's stated
  upper pressure, and a density above $N(p_\text{max}, T)$ is refused as outside it.

A bracketed solve (Brent's method, or Newton safeguarded by bisection) on the forward law
itself converges on this bracket from any start. The start is the backend's own
(density, temperature) evaluation. CoolProp's molar density is converted with each
fluid's own gas constant, $\rho_m = (N/N_A)(R_\text{CODATA}/R_\text{EOS})$, so that the
inverse keeps optmat's convention, $N = p/(Z k_B T)$ with CODATA's $k_B$
(`optmat.gas.eos`), rather than the equation's mass density. The start is then within the
EOS's own iteration tolerance of the root, and the solve takes one or two steps. The
fixed point $p_{k+1} = N k_B T Z(p_k, T)$ is not used: it converges linearly, with ratio
$|\partial\ln Z/\partial\ln p|$ at the root, which is not small near saturation.

The tolerance is set on the round trip, $|N(p(N,T),T)/N - 1|$, below a stated bound
(§8); a solve that does not reach it raises `StateError` rather than returning its last
iterate. Zero density, the vacuum limit, is solprop's to handle, and is refused
here.

#### 4.3 Transport of the dilute gases

The properties are:

- the dilute-gas shear viscosity $\eta(T)$, in Pa s: the zero-density limit of the
  viscosity, a function of temperature only, as kinetic theory and Wilke's rule assume;
- the ideal-gas molar isobaric heat capacity $c_p^\circ(T)$, in J mol⁻¹ K⁻¹: exactly the
  temperature derivative of the ideal-gas part of each gas's equation of state, not a
  value at a finite pressure. For the noble gases it is $\tfrac52 R$;
- the molar mass $M$, the one each gas's equation of state uses.

The viscosity at 0.1 MPa differs from the dilute limit by a few parts in 10⁴ (solprop
records $1.7612 \times 10^{-5}$ against $1.7599 \times 10^{-5}$ Pa s for N2 at 294 K,
1 bar and 1 mbar). That is below every source's stated uncertainty, so a table at
0.1 MPa may serve a dilute-gas value, with the offset stated in its record.

For a mixture of mole fractions $x_i$, Wilke's rule [2] gives

$$
\eta_\text{mix} = \sum_i \frac{x_i\,\eta_i}{\sum_j x_j\,\Phi_{ij}}, \qquad
\Phi_{ij} = \frac{\left[1 + (\eta_i/\eta_j)^{1/2}(M_j/M_i)^{1/4}\right]^2}{\left[8\,(1 + M_i/M_j)\right]^{1/2}}.
$$

The ideal mixture takes $c_{p,\text{mix}} = \sum_i x_i c_{p,i}$. Wilke's rule is a
kinetic-theory approximation. Its accuracy for the noble-gas and molecular mixtures of
interest is taken from the literature and stated on the manual page, and §5 lists the
reference.

What combines these properties for a flow problem stays in solprop: the dilute sound
speed $\sqrt{\gamma R T/M}$, the kinetic-theory mean free path and the Knudsen-regime
checks. solprop also offers a real-gas sound speed for a pure gas at a pressure, which it
takes from CoolProp directly. That is a property of the equation of state, so the
equation-of-state object of §6.5 gives it too, on the same phase domain; otherwise
solprop could not drop its direct pyfluids dependency, as SEP 3 intends.

#### 4.4 Air

For a named composition of mole fractions $x_s$, the air polarizability under solprop's
linear density scaling is the sum of its constituents:

$$
\alpha_\text{air}(\omega) = \sum_s x_s\,\alpha_s(\omega), \qquad
\chi = N\,\alpha_\text{air}/\varepsilon_0.
$$

The recommended models of the constituents give $\alpha_s$, so the air model is
consistent with the populations a propagation code carries. It has a derivative by
§4.1, and needs no second source.

**Which composition.** An air model is defined for a named composition, never for "air"
in general. Three are named, all from CIPM-2007 [1]:

- `"dry air (N2, O2, Ar)"`: N2, O2 and Ar renormalised to one (99.957 % of the
  molecules). It is the composition `raman_response("dry air")` already uses, and the
  default air model;
- `"dry air (N2, O2)"`: N2 and O2 renormalised, matching solprop's argon-free recipe;
- `"dry air (CIPM-2007)"`: all twelve constituents as printed. Each has an optmat
  model, so this composition needs no omission; it is renormalised by its printed sum,
  $1 - 2.9 \times 10^{-6}$, and says so.

Each named composition resolves to one model with its own identity, and a Raman response
and an air model built on the same name use the same fractions.

A refractometric air formula is a different, independent source. That is Ciddor's
(curated in part [3]), or Edlén's or Börzsönyi et al.'s, both curated. It therefore
serves as the validator. The comparison keeps two things apart:

- **omission**, the difference between a reduced composition and the full one,
  computed within optmat exactly. At `ea0dcaf` the N2/O2/Ar air's polarizability is
  1.9 × 10⁻⁴ below the full composition's from 400 to 1600 nm, and the N2/O2 air's
  1.5–1.7 × 10⁻⁴ above it;
- **model disagreement**, between the full composition and the formula at the
  formula's own CO2 fraction and state. It includes the refractivity virial of the
  mixture, which optmat does not model (solprop's SEP 11), and is bounded in §8.

**The validity of a composition.** An air model's qualified intervals are the
intersection of its constituents' intervals.
- Today every constituent has one interval, so this is simply the narrowest range. N2,
  O2 and Ar have no bands, so the two reduced compositions keep their present range for
  good.
- After OEP 3, the full composition's trace molecules bring their gaps. The
  intersection would then warn across CO2's 4.3 µm band at 400 µmol/mol, which is right
  in kind but too wide: a band's absorption scales with its gas's mole fraction.
- OEP 3 states when a mixture may widen the intersection (its §13 Q13). All three of
  these must pass at the mixture's mole fractions: the absorption weighted by mole
  fraction, no constituent pole in the wider interval, and each constituent's
  representation bound times its mole fraction within the tolerance. It is OEP 3's to
  supply, for air and for OEP 4's mixtures alike. Until then the intersection is the
  conservative rule.

**The density of air.** The air model's number density comes from the dry-air
equation of state CoolProp carries as a pseudo-pure fluid, which includes the
unlike-pair interactions. The sum of the constituents' densities at their partial
pressures, solprop's mixture rule, is lower by $1.9 \times 10^{-4}$ at 1 bar,
$1.8 \times 10^{-3}$ at 10 bar and $6.0 \times 10^{-3}$ at 50 bar (294 K, N2/O2/Ar,
against the pseudo-pure fluid's $Z$). The pseudo-pure fluid's own composition is to be
taken from its paper and stated. solprop's `Gas.Air` uses the pseudo-pure equation, and
SEP 3 states its difference from the partial-pressure rule of `Mixture.air` (§13, Q29).

#### 4.5 The Raman response and its groups

In optmat's convention ({doc}`../../raman_response`), following Hellwarth, the delayed
nonlinear polarisation is

$$
P_i(t) = \varepsilon_0 N \int_0^\infty R_{ijkl}(\tau)\, E_j(t)\, E_k(t-\tau)\, E_l(t-\tau)\, d\tau,
\qquad
R_{ijkl} = a(\tau)\,\delta_{ij}\delta_{kl} + \tfrac12\, b(\tau)\,(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}).
$$

A line's isotropic share enters $a$ only. Its anisotropic share $h_2$ enters with
$a = -h_2/2$ and $b = 3h_2/2$.

A propagation code working in a trace and traceless basis needs this written out for each
part separately:

- **The isotropic part:** $R^{(0)}_{ijkl} = h_0\,\delta_{ij}\delta_{kl}$.
- **The anisotropic part:**
  $R^{(2)}_{ijkl} = h_2\left[\tfrac34(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}) - \tfrac12\delta_{ij}\delta_{kl}\right]$.

Contracting with a linearly polarised field ($i=j=k=l$) gives $h_0$ and $h_2$ with unit
weight, as the scalar response $h = a + b$ requires. That check alone cannot tell the two
parts apart: an anisotropic line replaced by an isotropic one of the same amplitude
leaves it unchanged. The cross-polarised element does: a probe along $x$ after a pump
along $y$ sees $R_{xxyy} = a$, which is $h_0$ for the isotropic part and $-\tfrac12 h_2$
for the anisotropic part.

solprop's Raman leaf weights each term in the transverse trace and traceless basis:
$(1, 0)$ for an isotropic row and $(\tfrac14, \tfrac34)$ for a rank-two row. That is this
decomposition written for a transverse drive $Q$ embedded in three dimensions:
$\tfrac12 h_2(3Q - \operatorname{tr}Q\,I)$ has trace part $\tfrac14 h_2 \operatorname{tr}Q$
and traceless part $\tfrac32 h_2 (Q - \tfrac12 \operatorname{tr}Q\,I)$, which gives
$h_2(\tfrac14 + \tfrac34) = h_2$ copolarised and $h_2(\tfrac14 - \tfrac34) = -\tfrac12 h_2$
cross-polarised. The manual page derives the mapping, and §8 tests both elements.

Each line $k$ contributes a damped oscillator,
$A_k \sin(\Omega_k(N)\,t)\, e^{-\Gamma_k(N)\,t}$. With $N$ the total number density in
m⁻³ and $N_\text{amg}$ the amagat, $2.686780 \times 10^{25}$ m⁻³
(`optmat.gas.linewidths`):

- the decay rate is $\Gamma_k = \pi\,\Delta\nu_k$, where the full width at half maximum
  in Hz is $\Delta\nu_k(N) = c_k + \min(d_k N_\text{amg}/N,\ \Delta\nu_{D,k})
  + \beta_k N/N_\text{amg}$:
  a constant, a Dicke term capped at the forward-scattering Doppler width, and
  broadening (TC-687, TC-688). A line without a Dicke term has $d_k = 0$ and no cap;
- the angular frequency is $\Omega_k(N) = \omega_k + 2\pi s_k N/N_\text{amg}$, with
  the shift $s_k$ in Hz/amagat.

A **damping group** $g$ is a set of lines that share the width law and the shift at every
density: $(c, d, \beta, \Delta\nu_D, s)$ are equal for $k \in g$. That is how
`raman_response` numbers the groups, so the equality is exact. The group's response is

$$
R_g(t; N) = e^{-\Gamma_g(N)\,t} \sum_{k \in g} A_k \sin\!\big((\omega_k + 2\pi s_g N/N_\text{amg})\,t\big).
$$

The sum over groups equals the full response exactly. Where $s_g = 0$ the row of
amplitudes and frequencies is fixed and only the decay path depends on density: that is
the exact property solprop's term mapping relies on (D2). Where $s_g \ne 0$, the row
itself depends on density, and the mapping is exact only at a fixed density. At `ea0dcaf`
one group has a shift: H2's Q lines, at 98.4 MHz/amagat at 295 K. That is the
restriction SEP 3 places on H2 along a density gradient.

**The small-density limit.** Every term of the width law is finite for $N > 0$, and the
law has a finite limit as $N \to 0$: the capped Dicke term tends to its cap, so
$\Gamma_g(0^+) = \pi(c_g + \Delta\nu_{D,g})$ for a group with a Dicke term, and
$\pi c_g$ otherwise. For the curated gases at 295 K that limit is zero for every line of
N2, O2, dry air, N2O and CO2, and for some of D2's (widths of pure broadening), and
finite for H2 (at least 153 MHz), CH4 and SF6. There is therefore no lowest density at
which the law stops being defined, and no `minimum_density`. Three things are kept apart:

- **the limit**, which the group view gives at $N = 0$ exactly;
- **the validated range**, `density_range` (TC-688): below its lower end the lines are
  no longer independent collisional Lorentzians, and `fwhm`, `decay_rates` and
  `response` warn, as now;
- **a numerical clamp**, if a consumer wants one, which is the consumer's.

A response built by `from_lines` (§6.3) with a Dicke term and no cap would have an
infinite limit; the constructor refuses it.

#### 4.6 The energy per ionisation event

A rate model sums channels $c$, each with ionisation energy $I_c$ and rate
$\Gamma_c(E)$ at field amplitude $E$, the channel's weight included. The energy removed
from the field per ionisation event, consistent with the rate, is the rate-weighted mean

$$
\bar I(E) = \frac{\sum_c I_c\,\Gamma_c(E)}{\sum_c \Gamma_c(E)} = \frac{L(E)}{\Gamma(E)},
\qquad L(E) = \sum_c I_c\,\Gamma_c(E).
$$

$L$ is the **energy-loss rate** per particle, in W. It is the primary quantity, for three
reasons:

- **Averaging.** A consumer that averages the rate over a cycle or along a polarisation
  ellipse must average $L$ the same way. The averaged loss per event is
  $\langle L\rangle/\langle\Gamma\rangle$, which is not $\bar I$ at the peak field once
  $\bar I$ depends on the field. A consumer can integrate $L$ in the same cycle moment
  as the rate; one whose kernel carries a single energy per population (solprop's does)
  instead uses a constant energy, and $L$ and the channel rates let it state that
  approximation's error over the fields of its run.
- **Zero field.** $L$ and $\Gamma$ are both exactly zero at zero field, where $\bar I$ is
  0/0. $\bar I$ is defined there by its limit, the lowest channel's energy, since the
  tunnelling exponent $-2(2I_c)^{3/2}/(3|F|)$ suppresses the higher channels fastest.
  It is computed from the logarithms of the channel terms, so that it stays defined
  where every rate underflows.
- **The averaging factor.** The atomic default's "consistent" averaging multiplies the
  channel sum by one factor at each field, computed from the sum's slope and curvature
  in $\ln F$. The channel rates $\Gamma_c$ are the channel terms of that sum, each
  multiplied by the same factor, so $\sum_c\Gamma_c = \Gamma$ exactly and the factor
  cancels in $\bar I$. Separate single-channel models, each with its own factor, sum to
  within 0.1 % of the true rate for Ar at 800 nm, and change $\bar I$ by at most 2 meV.

For one channel, $\bar I$ is constant: for the molecular defaults it is the energy the
rate's exponent uses, vertical for every default at `ea0dcaf`. The alternative, the
adiabatic (ground-state) energy, is the thermodynamic threshold, but it disagrees with
the rate's own exponent (§7.8). For the atoms, $\bar I$ moves with the field between the
fine-structure channels (Ar at 800 nm: 15.797 eV at 5 × 10⁹ V/m, 15.814 eV at
5 × 10¹⁰ V/m, against channels of 15.760 and 15.937 eV). The movement is 0.25 % for Ar,
and up to about 5 % for Xe, whose channels are 1.31 eV apart.

The rates are for linear polarisation. The instantaneous rate along an elliptical field
trajectory is a quasi-static approximation, and is not validated by optmat.

#### 4.7 The silica nonlinear bundle

In the $n_0 = 1$ convention of a propagation code's nonlinear Schrödinger equation, the
nonlinear response is

$$
R(t) = (1 - f_R)\,\delta(t) + f_R\, h_R(t), \qquad \int_0^\infty h_R\,dt = 1,
$$

with the total $n_2$ from the Kerr fit. The consistent triple is $(n_2, f_R, h_R)$.

- **Paired:** the fraction is the Kerr fit's own. For Hollenbeck and Cantrell's shape,
  $f_R = 0.219$ (TC-696).
- **Unpaired:** the shape carries its paper's fraction: 0.245 for Lin and Agrawal's,
  0.222 for Blow and Wood's.

optmat's `raman_coefficient` and `electronic_susceptibility` include $n_0^2$. A
propagation code in the $n_0 = 1$ convention takes $n_2$, $f_R$ and $h_R$, not those.

**Which part a bundle holds.** A triple that is not paired cannot keep both the Kerr
fit's total $n_2$ and its electronic part $(1 - f_R^\text{fit})n_2$. Every bundle holds
the **total** long-pulse $n_2(\lambda)$ of `kerr/SiO2`, the quantity the $n_2$ data
measure and the one a GNLSE user supplies (SEP 3, D4), and takes its shape's own $f_R$.
Its electronic part is then $(1 - f_R)n_2$. The alternative, holding the electronic part,
gives a total of $n_2^{(e)}/(1 - f_R)$, 1.035 $n_2$ for Lin and Agrawal's fraction: that
is TC-696's "1.03 $n_2$".

**The gain check, run** (TC-734). At the Kerr fit's $n_2(526\text{ nm}) = 2.79 \times 10^{-20}$
m²/W, against Stolen and Ippen's $1.5 \times 10^{-11}$ cm/W ± 10 % at 330 cm⁻¹
(TC-693):

| Bundle | $f_R$ | Gain at 330 cm⁻¹ (10⁻¹¹ cm/W) | Against the measurement |
|---|---:|---:|---:|
| Hollenbeck, paired (the scalar default) | 0.219 | 1.35 | −9.7 % |
| Lin–Agrawal, scaled to the fit's $f_R$ | 0.219 | 1.23 | −17.8 % |
| Lin–Agrawal, its own $f_R$, total $n_2$ held | 0.245 | 1.38 | −7.9 % |
| Lin–Agrawal, its own $f_R$, electronic part held | 0.245 | 1.43 | −4.7 % |

Lin and Agrawal's normalised shape has an imaginary part 9 % smaller than Hollenbeck
and Cantrell's at 330 cm⁻¹, and a lower peak (1.71 against 1.91 × 10⁻¹¹ cm/W at the fit's
$f_R$, both near 435 cm⁻¹), so the same $f_R$ gives 9 % less gain. The fit's $f_R$ was set by this gain through Hollenbeck's shape (TC-696), so the
check is a consistency check for another shape, not an independent validation; the
independent validator of $f_R$ is Hellwarth et al.'s nuclear share, 0.21 ± 0.03
(TC-694). Lin and Agrawal's own 0.245 lies 1.0 standard deviation from the fit's
0.219 ± 0.026. §7.10 draws the consequence.

#### 4.8 Identity

The identity of a model is a digest that changes whenever the model's numbers can
change. The converse is not promised: an edit that leaves every number unchanged may
still change it. It is the SHA-256 [5] of a canonical record containing:

- the class's qualified name;
- the `EVALUATION_VERSION` integers of the class and of every shared evaluator module it
  evaluates through (`optmat.gas._tunnelling` for the rates, for example), each raised
  whenever that code changes its output;
- every constructor parameter, including the inputs that fix a state: temperature,
  wavelength, composition, population floor;
- **the values evaluation reads**, not the files that produced them: a fit's stored
  terms and offsets rather than its configuration's bytes, a measured spectrum's grid
  and values, a response's line arrays;
- **the declared validity**: the qualified intervals (§4.1), which decide where a model
  warns, and any stated error a model returns beyond its covariance (OEP 3's
  representation bound on the infrared terms);
- **the dependencies**, by their own digests (below).

**Dependencies.** A model's numbers may depend on another component. Two rules make the
record complete:

- **Resolve at construction.** A value obtained from another component and constant
  for the object's life is resolved once, at construction, and enters the record as a
  value: a published formula's reference density and refractivity-virial ratio, for
  example, which today are looked up (and cached) at each evaluation.
- **Bind what varies.** A component consulted at evaluation with a varying argument is
  bound to the object at construction, and its digest enters the record: the silica
  response's `raman_coefficient` reads `fit/SiO2`'s index and `kerr/SiO2`'s $n_2$ at the
  wavelength asked for; the bundle holds its Kerr model and response; an air model holds
  its constituents' models.

A backend computed at run time is a dependency like any other. The equation of state's
identity records CoolProp's version, the fluid, and the fluid's reference keys.

**Immutability.** A frozen dataclass can still hold writeable arrays and dicts, and the
registries return shared cached objects. Every model of the surface therefore owns its
storage: arrays are copied at construction and marked read-only
(`flags.writeable = False`), and mappings are copied into read-only mappings. A caller
that keeps and changes the arrays it passed changes nothing, and an attempt to change the
object's raises. The digest is then computed once and cached, and stays true for the
object's life.

**What the tests can prove.** Golden values, with the digest recorded beside them, fail
when a class's output changes at the sampled points without its `EVALUATION_VERSION`
being raised. That is a tripwire, not a proof: an output change away from the samples can
escape it. A consumer that cannot accept that risk also compares the version (§7.2), as
solprop SEP 3's first cut does.

The record is serialised as JSON:

- keys sorted, with no insignificant whitespace;
- UTF-8;
- floats written with Python's shortest round-trip `repr`, which is identical across
  IEEE-754 platforms;
- non-finite floats written as the strings `"inf"`, `"-inf"` and `"nan"`;
- arrays as lists;
- enums by value.

The package version is recorded beside the digest, not in it (§7.2). It is the full
version hatch-vcs gives, `importlib.metadata.version("optmat")`: a release number for a
tagged build, and the release, distance and commit (such as `0.0.1.dev264+g529dda475`)
otherwise. For an editable install that string is fixed when the package was installed,
not when it was last edited: the working tree of this review, at `ea0dcaf`, reported
`g529dda475`. The version is therefore informational, and the digest is the identity.

#### 4.9 The reference temperature of a polarizability

Every polarizability model and the refractivity virial are independent of temperature:
a state's temperature reaches the index through the density and Lorentz–Lorenz only
({doc}`../../conventions`, *Temperature*). The molecular fits declare the temperature
their polarizability refers to, 293.15 K for N2, O2, CO2, N2O, H2 and D2, and data
measured elsewhere are brought to it where a coefficient is measured (TC-730). A model
evaluated at another temperature does not follow it. The measured coefficients are, per
K: 1.18(4) ppm for N2 (`egan_2023`); 2.11 ppm for O2 and 4.22 ppm for CO2 at 633 nm,
each to about 9 % (`egan_2026`). H2's rotational average moves 0.03 % between 273 and
296 K. At 20 K from the reference the omission is about 24 ppm (N2), 42 ppm (O2) and
85 ppm (CO2) of $n - 1$.

The surface exposes, on each gas dispersion model, its reference temperature where its
fit declares one (none for the atoms, whose polarizability has no thermal average) and,
where one is curated, the
measured coefficient with its uncertainty, wavelength and source. Both enter the
identity's record. Evaluation is unchanged (§7.13).

The molecular fits' fixed infrared terms have a temperature of their own: HITRAN's
intensities refer to 296 K, and they include the lower-state populations. The metadata
records it beside the polarizability's reference temperature, since they differ by 3 K.
It matters most in the infrared windows of OEP 3, where hot bands (CO2's laser bands
among them) change strongly with temperature.

### 5. Primary sources and data

| Needed for | Source | State |
|---|---|---|
| Dry-air compositions | Picard et al., CIPM-2007 [1] | Curated (`picard_2008`) |
| Air validator | Ciddor [3]; Edlén; Börzsönyi et al. | Ciddor 2002 curated. Ciddor's 1996 equations [6] not yet curated. Edlén 1966 and Börzsönyi 2008 curated |
| Density of air | CoolProp's pseudo-pure dry air [9], Lemmon et al.'s equation [14] (CoolProp's reference key `Lemmon-JPCRD-2000`) | Curated (`lemmon_2000`, TC-741): Jones's composition, N2 0.7812, O2 0.2096, Ar 0.0092; 0.1 % in density (k = 2) to 873 K and 70 MPa |
| Viscosity of Ne, Kr, Xe, N2O | NIST Chemistry WebBook, SRD 69 [7], tables at 0.1 MPa, which evaluate REFPROP 10's extended-corresponding-states models; the models are documented by Huber [8] | Huber's report is in hand (public). It gives, for each gas, the equation of state, the data fitted and the uncertainty: Ne (§6.9.1: 2 % in the supercritical gas at atmospheric pressure), Kr (§6.9.2: 3 % in the gas at atmospheric pressure), Xe (§6.9.3: 0.7 % in the vapour), N2O (§6.11.2: 5 % in the gas; solprop's table records the WebBook's 10 %). The WebBook tables need their own record and extraction |
| Viscosity of D2 | Not yet identified | Huber's report has no deuterium model. solprop's table records that REFPROP evaluates the normal-hydrogen correlation of Muzny, Huber and Kazakov (2013) for deuterium. A primary for D2 is sought in Phase 2c (§13, Q14) |
| Viscosity of the other gases | CoolProp [9] (through pyfluids), and the primary correlation it implements for each | Not yet curated. CoolProp 8.0.0 names each correlation by a reference key (N2: `Lemmon-IJT-2004`); the record names the correlation's paper per gas, not the library alone |
| Ideal-gas $c_p^\circ$ and $M$ of every gas | The ideal-gas part of each gas's equation of state in CoolProp [9] | Not yet curated. CoolProp names each equation by a reference key (Ne: `Thol-JPCRD-2019-Neon`; Kr, Xe, N2O: `Lemmon-JCED-2006`; D2: `Richardson-JPCRD-2013`); the record names the paper per gas |
| Mixture viscosity | Wilke [2] | Curated (`wilke_1950`, TC-737): §4.3's formula is his Eqs. (13)-(14), checked against the page and reproduced from his Table I, He–Ar, Ne–Ar and Ne–He among its pairs. Stated accuracy: 0.97 % for his binaries, 1.9 % for multicomponent mixtures |
| Silica gain check | Stolen and Ippen [4] | Curated (TC-693) |
| Silica shapes | Hollenbeck and Cantrell [10]; Lin and Agrawal [11]; Blow and Wood [12] | Curated |
| Tensor form | Hellwarth et al. [13] | Curated (TC-694) |
| Digest | FIPS 180-4 [5] | A standard; no record needed |

**Vocabulary.** Two new quantities, each documented on {doc}`../../data` before any data
are added (the data rules):

- the dilute-gas viscosity (Pa s, abscissa $T$);
- the ideal-gas molar heat capacity (J mol⁻¹ K⁻¹, abscissa $T$).

The molar mass is a gas property with its source.

**Validators.**
- Ciddor's air against §4.4's sum for the full CIPM-2007 composition.
- Hellwarth et al.'s nuclear share of fused silica's $n_2$, 0.21 ± 0.03 (TC-694), for
  each bundle's $f_R$. Stolen and Ippen's gain is an input to the Kerr fit (TC-696), so
  against it the bundles are checked for consistency, not validated (§4.7).
- For the viscosities, a second source where one exists, compared on the page. For Ne,
  Kr and Xe, Huber's report states that the models' dilute-gas values at 25 °C and
  0.1 MPa were adjusted to agree with the reference values of Assael et al. (2018, as
  the report cites them), which are therefore not independent of them.

**Uncertainty.** Each transport record carries its correlation's stated uncertainty. The
air model's uncertainty is the constituents' combined. None of this changes an existing
fit.

References marked † in §14 are not yet curated source records. Each is checked against
the paper before acceptance, and no value is taken from them until it is.

### 6. Proposed interface

Sketches only, fenced as text until they run.

#### 6.1 The public surface

Each module listed in §3 declares `__all__`, as do the modules the surface adds:
`optmat.gas.transport`, `optmat.gas.air`, and the core module holding `ModelIdentity`
and the digest. The surface covers, besides the functions, every type a consumer names
or receives:

- the model base classes, `optmat.gas.models.base.Model`, with its qualified
  `intervals` (§4.1), and the solid base class;
- the explicit constructors of §6.3 and the classes they are on (the rate models,
  `GasRamanResponse`, `SilicaRamanResponse`, the solid `FormulaModel` and
  `TabulatedModel`);
- the returned types: `Channel`, the Raman group record, `SilicaNonlinearBundle`, the
  equation-of-state and transport objects, `Constituent`, `ModelIdentity`;
- the exceptions and warnings of `optmat.exceptions`.

A consumer never needs a private import or a name the page does not list. A manual page,
*Using optmat from a propagation code* (§9), lists every supported name, its signature
and the meaning of every returned field. `tests/test_public_api.py` snapshots the signatures and dataclass
fields of that surface. A change to the surface fails until the snapshot is updated
deliberately in the same commit, with a changelog line. The policy: before 1.0, a minor
version may break the surface and says so in the changelog; a patch version never does.

#### 6.2 Identity and flags

```text
@dataclass(frozen=True)
class ModelIdentity:
    name: str            # "fit/Ar", "kerr/SiO2", "raman_response/N2@295K,800nm", "rate/default/Ar@800nm"
    kind: str            # "gas_refractive" | "gas_hyperpolarizability" | "solid" | "kerr"
                         # | "raman_response" | "ionisation_rate" | ...
    material: str        # "Ar", "SiO2", "dry air", or a custom species' name
    digest: str          # §4.8
    optmat_version: str  # informational; not in the digest
    sources: tuple[str, ...]
    status: str          # "validated" | "provisional" | "untested" | "comparison" | "user"
    status_reason: str
    modified: tuple[tuple[str, str], ...] = ()   # overrides applied (§6.3), sorted by name
    base: str = ""       # the digest of the curated model an override started from
    reference_temperature: float | None = None   # gas dispersion models (§4.9), K
    infrared_temperature: float | None = None    # the fixed infrared terms' (§4.9), K
```

**Which objects have an identity.** Every model, response, bundle and composition
object that a surface constructor returns has `.identity`: the gas dispersion,
hyperpolarizability, rate, Raman-response and air models; the solid, Kerr and silica
Raman models and the bundle; the equation-of-state and transport objects of §6.5. A
function of the surface that returns plain numbers (`number_density`, `pressure`,
`viscosity`, ...) delegates to one of those objects, and the page names which, so its
provenance is that object's identity.

The status values are:

- `validated`: the default, with a stated uncertainty;
- `provisional`: a known stopgap (N2O's Raman widths, the Q widths of CO2 and N2O, a
  D2 viscosity served by a surrogate);
- `untested`: no reference (SF6's rate);
- `comparison`: a literature form kept for comparison (`ppt_ilkov`);
- `user`: built from parameters a caller supplied (`from_parameters`, `from_lines`,
  `from_arrays`), or changed by an override. Its `status_reason` names the base and its
  status where there is one; a curated model's stated uncertainty does not carry over.

#### 6.3 Explicit parameters and overrides

```text
IonisationRateModel.from_parameters(
    species: str, wavelength, *, ionization_energy, orbital_l, electrons,
    structure_coefficient | structure="hartree", coulomb="static",
    static_correction="tong-lin", static_constants=(...), channels=None, citation)
MolecularRateModel.from_parameters(
    species, wavelength, *, form, ionization_energy, ionization_energy_kind,
    structure_coefficients | effective_charge, electrons, tong_lin_alpha, ..., citation)

model.with_parameters(**changes)          # any model of the surface

GasRamanResponse.with_parameters(width_scale=..., amplitude_scale=..., width_law={band: ...})
GasRamanResponse.from_lines(species, temperature, frequency, isotropic, anisotropic,
    width_constant, width_dicke, width_broadening, doppler_width, shift,
    damping_group, citation)

FormulaModel.sellmeier(material, coefficients, wavelength_range, reference_medium, citation)
TabulatedModel.from_arrays(material, wavelength, n, k, citation)

OscillatorModel.from_parameters(species, strengths, energies, *, offset=0.0,
    intervals, citation)          # a gas polarizability as a sum of poles, SI
```

`OscillatorModel.from_parameters` builds a gas dispersion model from printed
coefficients, such as a published table or OEP 2's compact models, with the qualified
intervals the source states. The gas constructors that take a range accept `intervals`
as well; a single `wavelength_range` is one interval.

`from_parameters` validates exactly as the curated path does, and the curated builders
are rewritten to go through it, so that the two cannot differ. Its status is `user`, and
its identity's material is the given name.

**Overrides rebuild; they do not replace.** `with_parameters` is not
`dataclasses.replace`. Each class declares its **inputs**, the parameters an override may
name, and `with_parameters(**changes)` calls the class's `from_parameters` with the
base's inputs updated by the changes. Everything derived from the inputs is therefore
recomputed, and every constraint between inputs is checked. The rules, class by class:

- **`IonisationRateModel`.** It stores nothing derived; κ, n*, the structure
  coefficient, the barrier-suppression field and the maximum field are recomputed at
  each evaluation. The inputs are the channels, `orbital_l`, `structure` and
  `structure_coefficient`, `electrons`, `coulomb`, `static_correction` with its
  constants, `channel_sum` and `averaging`. `ionization_energy=` is accepted as a
  convenience that moves every channel, keeping the fine-structure splittings. The
  constraints already checked stay (`structure` with its coefficient, a correction with
  its number of constants); the missing ones are added (orbital $l \ge 0$, energies and
  a Tong–Lin α above zero).
- **`MolecularRateModel`.** It stores two derived values, `structure_factors` (the
  orientation average of the structure coefficients) and `coefficient_squared` (the
  Hartree coefficient of the PPT and ADK forms, from the energy and the charge). They
  are not inputs, and an override naming them is refused. The structure enters as
  coefficients, and the average is rebuilt from them. The existing
  `with_ionization_energy` becomes an override through the same path.
- **`GasRamanResponse`.** The inputs of `raman_response` (species or composition,
  temperature, wavelength, population floor) are not overrides: a different state is a
  new call, which rebuilds the populations, strengths and widths. The overrides act on
  the assembled lines: `amplitude_scale`, and `width_scale`, which multiplies $c$, $d$,
  $\beta$ and the Doppler cap alike and so keeps the groups; and `width_law={band: (c, d,
  β, cap)}`, after which the damping groups and `density_range` are recomputed by the
  same rules as `raman_response` uses (the group key; TC-688's limits).
- **`SilicaRamanResponse`.** The shape's parameters and `raman_fraction`. Any change sets
  `paired` to False.
- **The solid models.** The formula's coefficients, the table's arrays and the range.
- **The gas dispersion models.** The oscillator strengths, energies and offset, or a
  formula's coefficients, and the intervals.
  - A molecular fit's fixed infrared terms are one input, `infrared`, replaced whole or
    scaled by `infrared_scale`. There are thousands of them (4296 for CO2 today, more
    with OEP 3's rules), and no one overrides them one at a time.
  - A fit's covariance does not carry over.

`at_wavelength` on the rate models is not an override: the wavelength is a state, as the
Raman response's temperature is. Every override records each changed input in
`identity.modified`, sets the status to `user` with the base's digest in `base`, and
drops the base's stated uncertainty (a rate's `uncertainty_factor` becomes NaN unless the
caller gives one). α and $N_\text{eff}$ of the atomic default were calibrated with its
structure coefficients, for example, so a changed coefficient leaves their uncertainty
unstated.

One-liners:

```text
rate = ionisation_rate_model("Ar", 800e-9).with_parameters(structure_coefficient=2.4)
n2 = raman_response("N2", 295.0, wavelength=800e-9).with_parameters(width_scale=0.8)
```

#### 6.4 The polarizability derivative

`Model.dpolarizability_domega(wavelength, *, warn=True)` on the gas model base class:
complex, SI (C m² V⁻¹ s), with the range semantics of `complex_polarizability` (§4.1).
Each class implements its own analytic form (§4.1); the base class has no numerical
fallback, so a new class without one fails its test rather than silently differencing.

#### 6.5 Equation of state and transport

```text
optmat.gas.eos.equation_of_state(gas) -> EquationOfState   # .identity: CoolProp version,
                                                           # fluid, reference keys
EquationOfState.number_density(state) -> float             # m^-3, the forward law
EquationOfState.pressure(number_density, temperature) -> float   # Pa, §4.2
EquationOfState.speed_of_sound(state) -> float             # m/s, real gas
EquationOfState.phase_domain(temperature) -> tuple[float, float]  # Pa, (0, p_sigma or p_max)
optmat.gas.eos.number_density(gas, state, *, ideal=False)  # as now; delegates
optmat.gas.eos.pressure(gas, number_density, temperature, *, ideal=False) -> float  # Pa

optmat.gas.transport.transport_model(gas) -> TransportModel   # .identity, .sources,
                                                              # .status, ranges by quantity
TransportModel.viscosity(temperature) -> float        # Pa s, dilute gas
TransportModel.heat_capacity(temperature) -> float    # J/(mol K), ideal gas
TransportModel.molar_mass -> float                    # kg/mol
TransportModel.temperature_range(quantity) -> tuple[float, float]   # K
optmat.gas.transport.viscosity(gas, temperature) -> float           # delegates
optmat.gas.transport.heat_capacity(gas, temperature) -> float       # delegates
optmat.gas.transport.molar_mass(gas) -> float                       # delegates
optmat.gas.transport.mixture_viscosity(composition, temperature) -> float   # Wilke
```

Every function above accepts an array of temperatures or densities where it makes sense,
and refuses what it cannot answer with a specific `OptmatError` (`StateError` outside a
range or the phase domain), as solprop refuses now.

**The support matrix.** Every one of solprop's twelve pure gases has an equation of
state, a heat capacity and a molar mass, from CoolProp's equations; the viscosity comes
from CoolProp for seven and from the WebBook tables for Ne, Kr, Xe and N2O, and D2's
source is Phase 2c's to settle (§5). The manual page tabulates, per gas and quantity, the
source, the range and the status. Air has an equation of state (§4.4) and, through
Wilke's rule, a viscosity; it has no transport record of its own.

#### 6.6 Compositions and air

The named, renormalised compositions `"dry air (N2, O2, Ar)"` and `"dry air (N2, O2)"`
come from CIPM-2007, with their omissions stated in the docstring, beside the full
`"dry air (CIPM-2007)"`. `composition("dry air")` keeps its present meaning, the full
printed table; the air model's default name is `"dry air (N2, O2, Ar)"` (§4.4).

```text
optmat.gas.air.air_model(composition="dry air (N2, O2, Ar)") -> AirModel
AirModel: the gas Model interface (complex_polarizability, dpolarizability_domega,
    intervals, identity), plus .composition (the Constituents) and .equation_of_state
```

An air model is registered once per named composition, so solprop's two recipes resolve
to two identities. `raman_response` accepts the two reduced names, `"dry air"` remaining
an alias of the first, and refuses the full composition, whose CO2 has no broadening
coefficients by N2 and O2.

#### 6.7 Raman groups

- `GasRamanResponse.line_id[k]`: a stable identifier such as `"N2/rotation/S/6"` or
  `"H2/v1/Q/1"`, unchanged across releases while the line exists.
- `GasRamanResponse.groups()`: for each damping group, in SI and with the conversions of
  §4.5 written once, in optmat:
  - the line indices, and their `line_id`s;
  - the shared decay-rate law, $\Gamma_g(N) = \gamma_c + \min(\gamma_d/N,\ \gamma_D)
    + \gamma_\beta N$ in 1/s with $N$ in m⁻³ (so $\gamma_c = \pi c$,
    $\gamma_d = \pi d N_\text{amg}$, $\gamma_D = \pi\Delta\nu_D$,
    $\gamma_\beta = \pi\beta/N_\text{amg}$);
  - the shared shift $\sigma_g = 2\pi s_g/N_\text{amg}$, in rad s⁻¹ m³, and a flag for $\sigma_g
    \ne 0$;
  - the zero-density frequencies $\omega_k$, rad/s, and the rows of isotropic and
    anisotropic amplitudes, m⁵ V⁻² s⁻¹;
  - `decay_rate(N)`, vectorised and defined for $N \ge 0$, with the limit of §4.5 at
    $N = 0$. It does not warn: the range warning is the response's, given once for the
    densities a consumer declares (`check_densities(N)`).
- No `minimum_density` (§4.5).
- Vectorising `response` and `spectrum` over an array of densities is desirable for
  solprop's SEP 10, but not required.

#### 6.8 Rate energies

On the atomic and the molecular rate models, vectorised over the field amplitude (V/m):

- `channel_rates(field)`: the instantaneous rate of each channel, shape
  `(n_channels, *field.shape)`, summing exactly to `instantaneous_rate(field)` (§4.6);
- `energy_loss_rate(field)`: $L = \sum_c I_c\Gamma_c$, J/s per particle, zero at zero
  field; the quantity a consumer averages beside the rate;
- `energy_per_event(field)`: $\bar I = L/\Gamma$, J, the lowest channel's energy where
  the rate is zero or underflows;
- `energy_kind`: `"vertical"`, `"adiabatic"` or `"not stated"`, the kind of the energies
  the exponent uses, as the curated dataset states it. A model does not call an energy
  adiabatic because it sits in that field: SF6's `adiabatic_ionization_energy` holds a
  threshold whose kind its source does not state, and CO's holds a vertical energy.

A single-channel model has one channel. The linear-polarisation statement of §4.6 goes
in the docstring of each rate method.

#### 6.9 The silica bundle

```text
silica_bundle(name="default") -> SilicaNonlinearBundle(
    kerr: KerrModel, raman_fraction: float, response: SilicaRamanResponse,
    paired: bool, identity: ModelIdentity, status: str)
```

The bundle is one object whose parts cannot be mismatched. It holds the Kerr fit's total
$n_2(\lambda)$ and its shape's $f_R$ (§4.7), and gives $n_2$, $(1 - f_R)n_2$, $f_R$ and
$h_R$. Its names:

- `"default"`: Hollenbeck with the fitted $f_R$, paired, copolarised only;
- `"tensor"`: the tensor default of §7.10;
- one bundle per silica response, under the response's name.

The docstring of `raman_coefficient` states the $n_0 = 1$ rule of §4.7.

#### 6.10 Solids

- Identity and status on every solid model.
- `out_of_range="hold" | "raise"` on the table models. The default is `"hold"`, as now,
  so a caller can demand refusal.
- A capability note: N-BK7's index is lossless, and its transmittance is curated
  separately.
- `BK7` and `N-BK7` stay distinct.

#### 6.11 Parity fixtures

Every `tests/data/solprop_*.json` is re-recorded against a named solprop main-branch
commit, recorded inside the fixture. Two fixtures are new:

- `solprop_transport.json`: viscosity, heat capacity and molar mass over each gas's
  range;
- `solprop_eos.json`: number density and the pressure round trip.

These are the baseline for solprop SEP 3's changed-results report. OEP 3 refits CO, CO2,
CH4 and N2O. Fixtures holding those gases are re-recorded again when its refits land,
and the report says so.

#### 6.12 Release

A first version is tagged once §12's Phases 1–4 have merged. A test builds the wheel and
the sdist, installs them into a clean environment without the paper library or network
access, and runs the quickstart and the §6.1 surface. That proves the data files ship.
Until the tag, solprop pins an optmat commit.

### 7. Implementation strategies

The choices below were weighed before this OEP was written. Each is a recommendation to
be settled in review (§13), not a decision.

#### 7.1 The public surface: a facade or per-module `__all__`

- (a) One facade module, `optmat.api`, re-exporting the surface.
- (b) Per-module `__all__`, plus the manual page and the snapshot test.

**Recommended: (b).** A facade duplicates names and drifts from the modules it wraps.

#### 7.2 What the digest covers

- (a) Parameters, data and a per-class evaluation version (§4.8).
- (b) Also the package version.

**Recommended: (a).** It is stable across unrelated releases, so a propagation code can
resume across an optmat release that did not touch the model. A consumer can still
refuse on a version change if it wants to (solprop SEP 3's first cut does). Golden-value
tests keep (a) honest where they sample: they fail when a class's output changes there
without its `EVALUATION_VERSION` being raised (§4.8). The version string of an editable
install is fixed at install time (§4.8), which is a further reason not to make the
identity depend on it.

`ModelIdentity` and the digest function belong in the core (`optmat`), since both media
need them. The core still imports no medium.

#### 7.3 The status vocabulary

**Recommended:** the five values of §6.2 plus a free-text `status_reason`.

- The Kr and Xe rates with the suggested barrier-suppression constant are `validated`:
  Wallace et al.'s yields test them, and their uncertainty factor of 2 says how well.
- SF6's molecular default is `untested`.
- `user` is separate from `comparison`: a literature form is optmat's curation of a
  published model, while a model built or changed by a caller is not curated at all.
  The first draft gave `from_parameters` the status `comparison`, which conflated the
  two.

#### 7.4 Custom species

- (a) Accept a plain string species name in `from_parameters` now.
- (b) Wait for the generalisation of the `Gas` typing (TC-317).

**Recommended: (a).** The identity's material is the given name, and the generalisation
is recorded as TC-317's follow-up, not blocking this.

#### 7.5 The override API

- (a) Keyword arguments named after the dataclass fields, validated by name.
- (b) Typed parameter objects.

**Recommended: (a).** They are the cheapest to use in one line, and the dataclass already
is the parameter object.

#### 7.6 The transport boundary

**Recommended:** optmat owns the material properties: viscosity, $c_p$, molar mass, their
ranges and the mixture rule. solprop keeps the formulas that combine them for its flow
problem (§4.3).

#### 7.7 The air model

- (a) The composition's summed constituent polarizabilities (§4.4).
- (b) Ciddor's or Edlén's formula as the propagation model.

**Recommended: (a),** for each named composition, with `"dry air (N2, O2, Ar)"` as the
default (§4.4). It is consistent with the constituents a propagation code carries, it
has a derivative and an identity, and it needs no second source. (b) is good for
refractometry, but it is a different source, and it would disagree with (a) at the
10⁻⁴ level in the refractivity, from the omitted constituents and the models alike.

#### 7.8 The energy per event

- (a) The rate-weighted channel energy (§4.6).
- (b) The adiabatic energy.

**Recommended: (a),** exposed as the energy-loss rate $L$ with $\bar I$ derived from it
(§4.6). It is the energy the ionising channels actually cost, and it is the same
threshold that enters the rate's exponent.

#### 7.9 The weak rate defaults

**Recommended:** keep the suggested α of Kr and Xe, and SF6's molecular prediction, as
the defaults, with their flags and uncertainty factors, so that a consumer can warn once.
SF6 is flagged `untested`.

#### 7.10 The silica tensor pairing

- (a) Keep `lin_2006` unpaired: Lin and Agrawal's shape with its own $f_R = 0.245$ and
  the Kerr fit's total $n_2$.
- (b) Add a paired variant, `lin_2006_paired`, at the fitted $f_R$, gain-checked against
  Stolen and Ippen (§4.7).
- (c) Refit the Kerr fit with the tensor shape in the gain's prediction: a second stored
  fit, whose $f_R$ is paired with Lin and Agrawal's shape.

The first draft recommended (b), with (a) as the fallback if (b)'s gain were "far" from
the measurement. The check has now been run (§4.7): (b) is 17.8 % below the measured
gain, 1.8 times its stated uncertainty, where the scalar default is 9.7 % below and
(a) 7.9 % below.

**Decided: (a)** as the tensor default (J. C. Travers, 2026-09-30; TC-734), with status
`validated` and the gain check recorded, and (b) not added. Its $f_R$ is 1.0 standard deviation from the fit's, its gain
is as close to the measurement as the scalar default's, and it holds the same total
$n_2$, so a scalar and a vector run agree in the long-pulse limit. They differ in the
split: the vector run's electronic part is 3.3 % smaller and its Raman part 12 % larger.
(c) would pair the shape exactly; by the ratio of the two shapes' responses at 330 cm⁻¹
its $f_R$ would come out near 0.24, so it would change little beyond a second stored fit
and its robustness runs. It stays available as a later refinement.

"Far" (the former Q15) is then settled by the numbers rather than a rule: one stated
uncertainty of the measurement would reject (b), and (a) is within it.

#### 7.11 Release

**Recommended:** version 0.1.0, tagged once Phases 1–4 have merged and the clean-install
test passes. Tagging, pushing and publishing are actions for the maintainers.

#### 7.12 Cost and dependencies

- **Cost.** Every new operation is cheap: a derivative, a few Newton steps, a table
  interpolation, a SHA-256 of small records. Digests are computed once per construction
  and cached on the frozen object.
- **Dependencies.** No new runtime dependency: the EOS already uses pyfluids, and
  `hashlib` and `json` are in the standard library.

#### 7.13 The reference temperature (the former Q28)

- (a) Expose each gas dispersion model's reference temperature, and its measured
  coefficient where one is curated, with uncertainty, wavelength and source, as metadata
  that enters the identity (§4.9). Evaluation is unchanged.
- (b) Also offer an opt-in correction linear in $T - T_\text{ref}$, refused where no
  coefficient is curated (CH4, CO, SF6) and needless for the atoms.
- (c) Apply that correction by default.

**Decided: (a) in this OEP; (b) deferred to its own change** (J. C. Travers,
2026-09-30), on the {doc}`roadmap <../../roadmap>`. (a) is what a
propagation code needs to see the convention and record it. (b) is not needed by
solprop's SEP 3, whose default temperature, 294 K, is 0.85 K from the reference (1 ppm
of $n - 1$ for N2). It also carries an approximation that needs its own validation: the
coefficients of O2 and CO2 are measured at 633 nm, and a linear law does not hold for
H2's rotational average far from room temperature, the case a cold Raman cell presents.
(c) would change recommended values and the solprop parity fixtures, and needs
coefficients for every molecule first; the α(ω, T) data of Hohm's dissertation are on
the roadmap's wanted list and not in hand. solprop's own reference-temperature rule (§2,
out of scope) is a different matter: it is solprop's convention for its derived
quantities, not the models' polarizability.

### 8. Verification and validation

| Item | Check | Tolerance set by |
|---|---|---|
| Surface (§6.1) | The snapshot test is in the gate. The page lists every name solprop's SEP 3 adapters use, and nothing private is needed to build any of them. Every gas model's `intervals` is its present range, one interval, and every stored model warns exactly where it did | Exact |
| Identity (§6.2) | Identical inputs give identical digests across processes (a subprocess test). Any parameter or data change changes the digest (a property test over the fields). Changing a dependency (the Kerr fit behind a silica response, a constituent behind an air model, the backend version behind an EOS) changes the dependent digest, and so does changing a model's intervals. Writing to an array or mapping of a model raises, and changing an array a caller passed changes nothing. Golden values, with their digests, fail when `EVALUATION_VERSION` is not raised | Exact |
| Overrides (§6.3) | `from_parameters` with the curated inputs reproduces the curated model exactly, for every gas, compared through the digest's parameter record. For every declared input of every class, `with_parameters(x=v)` equals `from_parameters` of the changed inputs, in its digest and its numbers, including the derived values (the molecular structure factors and Hartree coefficient; the Raman groups and `density_range` after a width law changes). A derived field is refused. An override sets `modified`, `base` and the status `user`. `from_lines` round-trips a curated response, and refuses a Dicke term without a cap | Exact, or ulp where arithmetic is reordered |
| Derivative (§6.4) | Every registry model agrees with a multi-step central difference of `complex_polarizability`, converging at second order, away from poles. O3 is tested inside its bands, between knots, outside the grid and at $\omega = 0$, in the real and the imaginary part; the knot convention is tested at a knot. He's retardation term is tested. The static limit is 0 | The difference's truncation and rounding error, stated in the test |
| EOS (§6.5) | The round trip $N \to p \to N$ holds on each gas's phase domain at 0.01 bar up to its saturation pressure (below $T_c$) or 100 bar (above), over its temperature range. The ideal limit is recovered as $p \to 0$. The saturated vapour's density, a density above it, a liquid state and a supercritical-liquid state are refused, forward and inverse. The inverse keeps the CODATA convention: its $N$ agrees with `number_density`, not with the backend's mass density | Solver tolerance, stated |
| Transport (§6.5) | Viscosities equal solprop's current values where the sources are the same, to their 1 bar–dilute offset. Differences are documented on the page where they are not. $c_p^\circ$ of the noble gases is $\tfrac52 R$. Wilke's rule reduces to the pure gas at $x_i = 1$, and is symmetric under relabelling. Out-of-range calls are refused | The sources' stated uncertainty |
| Air (§6.6) | Each air model equals its composition's weighted sum of the constituents' models, value and derivative. The two reduced compositions give solprop's two recipes to their rounding. The full composition agrees with Ciddor's formula within the constituents' combined uncertainty and the unmodelled mixture virial; the reduced compositions' difference from the full one is reported separately as omission. The air EOS is compared with the partial-pressure sum on the page. An air model's intervals are the intersection of its constituents' | Combined uncertainty; summation bound for the sums |
| Raman groups (§6.7) | Group members' width laws and shifts are equal (0 ulp) for every gas. The sum over groups equals `response` to summation roundoff at 100 densities. Line identifiers are unique and stable (snapshot). `decay_rate(0)` is the limit of §4.5, and `decay_rate(N)` tends to it as $N \to 0$. The copolarised element equals iso + aniso and the cross-polarised element iso − ½ aniso, line by line, and the $(1, 0)$ and $(\tfrac14, \tfrac34)$ weights reproduce both. The fixed-row property is asserted only for groups without a shift | Summation bound |
| Rate energies (§6.8) | `channel_rates` sum to `instantaneous_rate`. `energy_loss_rate` equals the channel sum of energy times rate, and is 0 at zero field. `energy_per_event` equals their ratio at every field, the lowest channel's energy at zero field and where the rates underflow, and for a single channel equals the rate's `ionization_energy` | Roundoff |
| Silica (§6.9) | Each bundle's total is the Kerr fit's $n_2$, and its parts add up to it. `paired` implies $f_R$ equal to the Kerr fit's. The gain check of §4.7 is recorded for every bundle | Stolen and Ippen's stated uncertainty |
| Solids (§6.10) | Identity tests. The raise option is tested | Exact |
| Release (§6.12) | The clean-install test | Pass/fail |

Every tolerance is justified in a comment in its test (`AGENTS.md`). Each phase ends with
its register entries, its changelog line and a green gate.

### 9. Documentation

- **A new page, *Using optmat from a propagation code*.** It covers:
  - the surface (§6.1);
  - the units: SI, vacuum wavelength, per particle;
  - which helpers use solprop's convention (`gamma3`) and which do not ($n_2$, $K$,
    $\chi_e$);
  - identities and the meaning of the digest;
  - qualified intervals, the envelope and the range warning;
  - the flags;
  - overrides and `from_parameters`, with the one-liners;
  - the Raman group view and the zero-density limit;
  - the rate energy per event;
  - the silica bundle.

  Its snippets are doctested, as the quickstart's are.
- **A transport section or page:** the sources, ranges and Wilke's rule, with the
  comparison against solprop's current values.
- **Updates to existing pages:**
  - {doc}`../../raman_response`: the tensor derivation of §4.5, groups, line identifiers,
    the zero-density limit;
  - {doc}`../../ionisation_rates` and {doc}`../../molecular_ionisation_rates`: channel rates,
    the energy-loss rate and energy per event, linear polarisation, flags;
  - {doc}`../../silica`: the bundles and their gain checks (§4.7);
  - {doc}`../../metals`: `out_of_range`;
  - {doc}`../../data`: the new quantities;
  - {doc}`../../conventions`: the reference temperature as exposed metadata (§4.9), and
    the equation of state's phase domain and inverse (§4.2).
- **Docstrings** with units and ranges on every new public name. Each page produced ends
  with a link to this OEP at Final.

### 10. Compatibility and impact

- **Numbers.** No recommended value of an existing model changes. These things do change:
  - SF6's rate becomes flagged `untested`;
  - the silica tensor default is named: `lin_2006` with its own $f_R$ and the Kerr fit's
    total $n_2$ (§7.10(a), TC-734). The scalar default is unchanged;
  - the forward density refuses a condensed state it now returns silently (§4.2). None
    of the 122 distinct gas states at which curated datasets are stated is condensed (a
    check run at `ea0dcaf`); a register entry records the finding and the change;
  - O3's first derivative (`dpolarizability_domega_au`, and `dispersion_au` at order 1)
    moves to the spline's own, by at most 4.6 × 10⁻⁵ of its maximum within O3's range
    (§4.1). No fit's objective uses it. The higher orders of `dispersion_au`, which
    `fit/O3`'s GVD uncertainty, its robustness record and its gallery page use, stay
    repeated gradients of the grid: a cubic spline's third derivative is piecewise
    constant, so they need their own treatment, recorded on the roadmap.
- **Interface.** The surface is declared and snapshotted, and future breaks are
  deliberate. `intervals` is new, with one interval per model, so no warning changes;
  `wavelength_range` keeps its value and is documented as the envelope (§4.1). Names become public by being added to `__all__`. No existing public name
  is removed.
- **Fits.** None is re-run by this OEP. §7.10(c), a second Kerr fit with Lin and
  Agrawal's shape, would be a separate decision.
- **Fixtures.** Re-recorded, with the solprop commit inside each.
- **Register.** Every source record, fit-adjacent decision and finding gets its entry,
  starting at the next free number at implementation.

### 11. Alternatives considered

- **A facade module** (§7.1). Rejected because it duplicates and drifts.
- **Versioned digests** (§7.2). They invalidate every downstream cache on every release.
- **Typed override objects** (§7.5). Heavier than keywords, for no gain.
- **Refractometric air** (§7.7). An inconsistent second source.
- **The adiabatic energy** (§7.8). Inconsistent with the rate's exponent.
- **Refitting the Kerr fit to a tensor shape** (§7.10(c)). A larger change than is
  needed.
- **Leaving transport in solprop.** It would keep sourced material data out of the
  register, contrary to SEP 3's aim.
- **Numerical derivatives everywhere** (§4.1). Cheaper to write, but less accurate, and
  unnecessary where the models are closed forms.

### 12. Phases

Each phase is independently mergeable with a green gate, its register entries and its
changelog line.

- **Phase 1 — Surface and identity** (§6.1, §6.2, §4.8). `__all__` and the snapshot
  test, `ModelIdentity` in the core with the digest and `EVALUATION_VERSION`, the flags,
  `Model.intervals` with one interval per model (§4.1), and golden-value tests. OEP 3's
  Phase 2 builds on it. Owned, read-only storage in every model of the surface, and
  the dependencies resolved or bound at construction. Every later phase needs it.
- **Phase 2 — The independent pieces**, in parallel after Phase 1, each on its own
  branch:
  - 2a: explicit parameters and overrides (§6.3);
  - 2b: the polarizability derivative (§6.4), and the reference temperature as metadata
    (§4.9);
  - 2c: the EOS object, its phase domain and inverse, and transport, with the new
    quantities and source records, including D2's viscosity (§4.2, §6.5);
  - 2d: compositions, the air models and the air EOS (§6.6);
  - 2e: Raman groups, line identifiers, the zero-density limit and the tensor derivation
    (§6.7);
  - 2f: channel rates, rate energies and flags (§6.8);
  - 2g: the silica bundles and the tensor default (§6.9, §7.10);
  - 2h: solids (§6.10).

  Phases 2e, 2f and 2g use 2a's `with_parameters` where they add overrides. They may
  merge before it and add the overrides afterwards.
- **Phase 3 — Parity fixtures** (§6.11), after Phase 2.
- **Phase 4 — The propagation-code page** (§9), after Phase 3.
- **Phase 5 — Release** (§6.12): the clean-install test. The tag itself is left to the
  maintainers.

### 13. Open questions

Each question below is kept with its outcome. **Resolved** means the text now answers it;
**Decided** records the maintainers' decision; **Open** is what remains before
acceptance; **Deferred** names the phase that settles it. The GPT-6 review of 2026-09-29
raised Q18–Q27; the review of 2026-09-30 (Opus 5.5) resolved them against the code at
`ea0dcaf`, and raised Q29. J. C. Travers accepted that review's recommendations on
2026-09-30.

**The recommendations of §7 (Q1–Q11).**

- **Q1** (§7.1, per-module `__all__`), **Q2** (§7.2, the digest without the version),
  **Q4** (§7.4, a plain species name), **Q5** (§7.5, keyword overrides), **Q6** (§7.6,
  the transport boundary), **Q7** (§7.7, the summed air model), **Q8** (§7.8, the
  rate-weighted energy), **Q9** (§7.9, the weak defaults kept and flagged) and **Q11**
  (§7.11, 0.1.0 after Phases 1–4). *Decided (2026-09-30):* confirmed as recommended;
  the review found no reason to change them. It sharpened how Q2, Q5, Q6, Q7 and Q8 are carried out: §4.8, §6.3,
  §4.3 with §6.5, §4.4 and §4.6. For Q4, the models' `gas` field, typed `Gas` today, widens to `Gas | str`.
- **Q3** (§7.3). *Resolved:* revised to five values, adding `user` for anything a caller
  built or changed.
- **Q10** (§7.10). *Decided (2026-09-30, TC-734):* the gain check has been run (§4.7),
  and the tensor default is (a), the unpaired `lin_2006` with the total $n_2$ held, in
  place of the paired variant. solprop SEP 3's D11 anticipated this fallback, and its
  text is brought into line.

**The first draft's questions (Q12–Q17).**

- **Q12.** Should acceptance wait for solprop SEP 3's? *Resolved: no.* Nothing here is a
  code path only solprop uses (§2), and every piece is useful to any propagation code.
  The order is: this OEP accepted; SEP 3 revised where this OEP's resolutions touch it
  (its §4.4 small-density clamp, its §4.6 energy loss, whose kernel carries one energy
  per population, D9's composition name,
  D11, the equation-of-state domain, and `Gas.Air`'s density, Q29); then SEP 3 accepted
  before its Phase 1.
- **Q13.** Does any gas model need a numerical derivative? *Resolved: none* (§4.1). Every
  registry model has an analytic form, O3's measured term included, branch by branch.
  `euv.ComplexModel` is outside the registry and the scope.
- **Q14.** Which correlation is behind each tabulated viscosity? *Resolved for Ne, Kr, Xe
  and N2O* (§5: Huber's report, sections 6.9.1–6.9.3 and 6.11.2), and for the gases
  CoolProp serves (each correlation's reference key). *Deferred to Phase 2c for D2:* no
  D2 model is in Huber's report, and solprop's table evaluates a normal-hydrogen
  correlation for it. Phase 2c seeks a primary for D2. If none is in hand, D2's viscosity
  is served from the WebBook table with the status `provisional` and the surrogate
  stated. The phase does not merge until one or the other is recorded in the register.
- **Q15.** How far may the paired gain be from the measurement? *Resolved* by running
  the check (§4.7, §7.10).
- **Q16.** Which version string? *Resolved:* the full hatch-vcs version, informational
  only, with the editable-install caveat (§4.8).
- **Q17.** The † references. *Resolved.* The bibliographic data of every †
  reference were checked on 2026-09-30 against its DOI record: [2], [5], [6], [8] and
  [9] through Crossref, and [7]'s DOI resolves to the WebBook. Huber's report [8] is
  public, and §5's statements from it were read in it. Nothing else is taken from a †
  reference before its source record exists. *Resolved (2026-09-30):* Wilke's paper [2]
  is curated (`wilke_1950`, TC-737). §4.3's formula is his Eqs. (13)-(14), checked
  against the page image, and it reproduces his Table I. Lemmon et al.'s air equation
  [14] is curated too (`lemmon_2000`, TC-741).

**The GPT-6 review (Q18–Q27).** *All resolved.*

- **Q18.** Identity dependencies and immutability: §4.8. Values resolved at
  construction, varying dependencies bound with their digests, the backend's version,
  owned read-only storage, and the promise weakened to "whenever".
- **Q19.** Overrides must rebuild derived values: §6.3. Declared inputs per class, a
  rebuild through `from_parameters`, derived fields refused, state changes left to the
  constructors, the status `user` with the base's digest, and the base's uncertainty
  dropped; §8 tests equivalence with a fresh build for every input.
- **Q20.** The small-density contract: §4.5 and §6.7. The width law has a finite limit
  at zero density, so `minimum_density` is withdrawn. The group view gives the limit,
  TC-688's range and its warnings are unchanged, and a clamp is the consumer's.
- **Q21.** Which air composition: §4.4 and §6.6. Three named compositions; the default
  is `"dry air (N2, O2, Ar)"`, the Raman response's; one identity per name; the Ciddor
  comparison separates omission from model disagreement.
- **Q22.** The derivative of the actual absorption evaluator: §4.1. The spline's
  derivative inside the grid, the differentiated trapezoid outside, the segment slope
  for the imaginary part with a stated knot convention; the existing
  `real_derivative_au` is not reused.
- **Q23.** The rate-energy domain and averaging: §4.6 and §6.8. The energy-loss rate
  $L$ is primary and is averaged like the rate; the channel rates sum exactly to the
  rate, the averaging factor being common to them; $\bar I$ takes the lowest channel's
  energy at zero field.
- **Q24.** Group units, shifts and the tensor: §4.5, §6.7 and §8. SI coefficients with
  the conversions written once; the fixed-row property only without a shift (H2's Q
  lines are the one shifted group); the cross-polarised element tested, and solprop's
  $(\tfrac14, \tfrac34)$ weights derived.
- **Q25.** The EOS inversion domain: §4.2. A declared phase domain shared by the forward
  law and its inverse, a bracketed solve on a monotone function, refusal of condensed
  states, and the CODATA convention kept.
- **Q26.** The silica fallback and the gain criterion: §4.7 and §7.10. Every bundle
  holds the total $n_2$; the check has been run; Stolen and Ippen's gain is a
  consistency check, and Hellwarth et al.'s share is the validator.
- **Q27.** The complete transport and surface contract: §4.3, §5, §6.1, §6.2 and §6.5.
  Heat capacities from each equation's ideal-gas part for all twelve gases; the
  transport and EOS objects carry the identities; the inventory lists the constructors,
  base classes and returned types. D2's viscosity remains, under Q14.

**Later questions.**

- **Q28.** The polarizability's own temperature dependence. *Decided (2026-09-30):* (a),
  the reference temperature and coefficients as metadata in the identity, now (§4.9);
  (b), the opt-in correction, as a later change with its own validation, on the roadmap
  (§7.13).
- **Q29.** The density of air. `Gas.Air` needs a number density. optmat's air model will
  carry CoolProp's pseudo-pure dry air, which includes the unlike-pair interactions;
  the partial-pressure sum solprop uses for mixtures is lower by 1.9 × 10⁻⁴ at 1 bar and
  1.8 × 10⁻³ at 10 bar (§4.4). *Decided (2026-09-30), in solprop's SEP 3:* `Gas.Air`
  uses the pseudo-pure equation, the better physics, and SEP 3 states its difference
  from `Mixture.air`'s partial-pressure rule. This OEP supplies both.
- **Q30.** What OEP 3's infrared windows need from this surface. *Resolved in the text
  (2026-09-30):*
  - validity as qualified `intervals`, one per model here, with `wavelength_range` as
    the envelope (§4.1, §6.1);
  - the intervals and any extra stated uncertainty in the identity (§4.8);
  - the infrared terms' 296 K beside the reference temperature (§4.9, §6.2);
  - an air model's intervals as the intersection of its constituents', until OEP 3
    supplies the checks under which a mixture may widen it (§4.4);
  - the infrared terms as one override input, and `OscillatorModel.from_parameters`
    with intervals (§6.3);
  - the molecular fixtures re-recorded after OEP 3's refits (§6.11).

  None of this changes a number or a warning now. The author asked for these additions;
  they are part of the human review of Part I.

### 14. References

1. A. Picard, R. S. Davis, M. Gläser and K. Fujii, "Revised formula for the density of
   moist air (CIPM-2007)," *Metrologia* **45**, 149 (2008). doi:10.1088/0026-1394/45/2/004
2. C. R. Wilke, "A viscosity equation for gas mixtures," *J. Chem. Phys.* **18**,
   517–519 (1950). doi:10.1063/1.1747673 (`wilke_1950`)
3. P. E. Ciddor, "Refractive index of air: 3. The roles of CO2, H2O, and refractivity
   virials," *Appl. Opt.* **41**, 2292 (2002). doi:10.1364/AO.41.002292
4. R. H. Stolen and E. P. Ippen, "Raman gain in glass optical waveguides," *Appl. Phys.
   Lett.* **22**, 276 (1973). doi:10.1063/1.1654637
5. † National Institute of Standards and Technology, *Secure Hash Standard*, FIPS
   PUB 180-4 (2015). doi:10.6028/NIST.FIPS.180-4
6. † P. E. Ciddor, "Refractive index of air: new equations for the visible and near
   infrared," *Appl. Opt.* **35**, 1566 (1996). doi:10.1364/AO.35.001566
7. † P. J. Linstrom and W. G. Mallard (eds.), *NIST Chemistry WebBook*, NIST Standard
   Reference Database 69, thermophysical properties of fluid systems. doi:10.18434/T4D303
8. † M. L. Huber, *Models for viscosity, thermal conductivity, and surface tension of
   selected pure fluids as implemented in REFPROP v10.0*, NISTIR 8209 (NIST, 2018).
   doi:10.6028/NIST.IR.8209. Public, and read for §5.
9. † I. H. Bell, J. Wronski, S. Quoilin and V. Lemort, "Pure and pseudo-pure fluid
   thermophysical property evaluation and the open-source thermophysical property
   library CoolProp," *Ind. Eng. Chem. Res.* **53**, 2498–2508 (2014).
   doi:10.1021/ie4033999
10. D. Hollenbeck and C. D. Cantrell, "Multiple-vibrational-mode model for fiber-optic
    Raman gain spectrum and response function," *J. Opt. Soc. Am. B* **19**, 2886
    (2002). doi:10.1364/josab.19.002886
11. Q. Lin and G. P. Agrawal, "Raman response function for silica fibers," *Opt. Lett.*
    **31**, 3086 (2006). doi:10.1364/ol.31.003086
12. K. J. Blow and D. Wood, "Theoretical description of transient stimulated Raman
    scattering in optical fibers," *IEEE J. Quantum Electron.* **25**, 2665 (1989).
    doi:10.1109/3.40655
13. R. Hellwarth, J. Cherlow and T.-T. Yang, "Origin and frequency dependence of
    nonlinear optical susceptibilities of glasses," *Phys. Rev. B* **11**, 964 (1975).
    doi:10.1103/PhysRevB.11.964
14. E. W. Lemmon, R. T. Jacobsen, S. G. Penoncello and D. G. Friend, "Thermodynamic
    properties of air and mixtures of nitrogen, argon, and oxygen from 60 to 2000 K at
    pressures to 2000 MPa," *J. Phys. Chem. Ref. Data* **29**, 331 (2000).
    doi:10.1063/1.1285884 (`lemmon_2000`)

[sep3]: https://github.com/LupoLab/solprop/blob/main/docs/dev/seps/sep-0003-optmat-materials.md

## Part II — Implementation record

### Phase 1 — Surface and identity

- **Scope:** Phase 1 of §12, 2026-09-30, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0001-identity`, merged into `main`.
- **Delivered:**
  - `src/optmat/identity.py` (new, core): `ModelIdentity`; the canonical record and its
    SHA-256 (§4.8); `Identified`, the base of every model of the surface, with its
    status fields, `IDENTITY_KIND`, `EVALUATION_VERSION`, `EVALUATES_THROUGH`, the
    cached `identity` and `identity_record()`; `FrozenMapping` and `frozen`, the owned,
    read-only storage; `parse_status` and the status vocabulary; `source_ids`.
  - `src/optmat/_validity.py` (new, core, private): the qualified intervals (§4.1).
  - `intervals` on the gas and solid model bases, one interval per model; the range
    warning goes through it.
  - `Identified` on `optmat.gas.models.base.Model` and its subclasses,
    `MeasuredAbsorption`, `HyperpolarizabilityModel`, `IonisationRateModel`,
    `MolecularRateModel`, `GasRamanResponse` (with its new `material` and
    `population_floor` fields), the solid models, `KerrModel` and
    `SilicaRamanResponse`. A module-level `EVALUATION_VERSION` in the six evaluator
    modules they name.
  - The dependencies resolved or bound once: the virial ratio of every gas dispersion
    model, the reference density of a formula (`reference_density`), silica's
    `index_model` and `kerr_model`, a solid formula's `reference_model`.
  - The statuses: `gas/models/_defaults.py` and `solid/models/_defaults.py` (the
    default rule, moved out of the registries), the rate and Raman registries, and
    `status` keys in `data/gas/molecular_ionisation_rates.toml` (SF6) and
    `data/gas/raman_response.toml` (CO2, N2O).
  - Tests: `tests/test_identity.py`, `tests/test_model_identity.py`,
    `tests/test_golden_models.py` with `tools/golden_models.py` and
    `tests/data/golden_models.json`, `tests/test_public_api.py` with
    `tools/public_api.py` and `tests/data/public_api.json`.
  - The manual page {doc}`../../identity`, and the core's entry in the API reference.
- **Deviations from the proposal:**
  - The values resolved from another component are resolved once per object on first
    use, not at construction (§4.8). Eager resolution would make building the gas
    registry import CoolProp and fail without the optional pyfluids, and would make
    every solid model load the gas fits. The value is still fixed for the object's life
    and in its record (TC-745).
  - §4.8's example of a shared evaluator, `optmat.gas._tunnelling`, is not one: it
    computes stored structure factors at construction. The rate models name
    `ionisation_rates`, where the atomic evaluator lives (TC-745).
  - The registry functions still return their shared, cached dictionaries. The models
    in them are read-only; a caller can still add to or remove from a dictionary.
  - `__all__` needed no change: every module of §3 already declared one. The snapshot
    covers the 27 modules of `tools/public_api.py`.
  - The page *Using optmat from a propagation code* is Phase 4's; this phase adds its
    first part as {doc}`../../identity`.
- **Decisions:** the status of every curated model (TC-746). `Model.kind` keeps its
  meaning, the provenance, and `ModelIdentity.kind` names the quantity; the page says
  so.
- **Register entries:** TC-745 to TC-746.
- **Verification:**
  - The full gate, `uv run poe check`, on the branch before the merge.
  - Golden values of 358 quantities were recorded on `main` before any change; after
    it every one is bitwise unchanged. Their digests are now recorded beside them, and
    the test fails when a value changes while its digest does not.
  - Digests are identical across processes and hash seeds (a subprocess test), and
    change with any field, an interval, an evaluator version or a bound component
    (§8's identity row). Every curated model's digest is distinct.
  - Storage is read-only after construction, `dataclasses.replace` and pickling, for
    every model of every registry.
  - Every model has one interval, its range; the warning's message for one interval is
    unchanged, and a NaN wavelength still does not warn.
  - Not run: `test-fits` and `test-library`, since no fit, extraction or curated value
    changed.
- **Sources checked:** none needed; no value was transcribed. The statuses cite the
  register entries that established them.
- **Limitations:**
  - The golden digests are recorded on one platform. A platform whose arithmetic rounds
    a stored value differently gives a different digest, as it should, and would fail
    the golden test.
  - The flags a consumer reads are the statuses; the SF6 rate and the CO2 and N2O
    Raman responses are the only non-validated defaults.
- **Open questions:**
  - Q31. Should the registry functions return read-only mappings? It would change their
    return types. For Phase 2a, with the constructors.
  - Q32. Where the golden digests should be recorded once CI runs on other platforms.
    For Phase 5.
- **Handoff:** Phase 2 builds on `Identified`. A new model class declares
  `IDENTITY_KIND`, its own `EVALUATION_VERSION`, `_identity_label()` and, where it
  reads another component, `_identity_extras()`, and calls `super().__post_init__()`.
  OEP 3's Phase 2 sets `intervals` on the molecular fits. Re-record
  `tests/data/golden_models.json` and `tests/data/public_api.json` with their tools
  when a model or the surface changes on purpose, with a changelog line.

### Phase 2b — The polarizability derivative and the reference temperature

- **Scope:** Phase 2b of §12, 2026-09-30, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0001-derivative`, merged into `main`.
- **Delivered:**
  - `Model.dpolarizability_domega`, with an abstract `_dcomplex_polarizability_domega`
    that every class implements. `OscillatorModel` (and so `FittedModel`), the gas
    `FormulaModel` and `DatabaseModel` do so analytically (§4.1).
  - `MeasuredAbsorption.complex_derivative_au` gives the derivative branch by branch.
  - `lorentz_lorenz.dpolarizability_dsusceptibility` gives the Jacobian of inverse Lorentz–Lorenz.
  - The reference temperature (§4.9): `FittedModel.reference_temperature`,
    `infrared_temperature` and `temperature_coefficient` (`TemperatureCoefficient`), in
    the identity's record and in `ModelIdentity.reference_temperature` and
    `infrared_temperature`; `data/gas/temperatures.toml`.
  - `tests/test_polarizability_derivative.py`; the derivatives added to the golden
    values; the models and conventions pages.
- **Deviations from the proposal:**
  - The measured coefficients are curated in a small table, `temperatures.toml`, beside
    the fits, rather than read from the fit configurations. Editing those would
    change the fits' digests and require refits. A test checks the two agree.
  - The coefficient's wavelength comes from its source record (632.9908 nm for all three).
- **Decisions:** O3's first derivative moves to the spline's own, as §10 states (TC-757).
- **Register entries:** TC-757.
- **Verification:**
  - The full gate on the branch.
  - Every registry model against a four-point difference, and O3 branch by branch, as TC-757 records.
  - The golden values are unchanged; their digests moved with `OscillatorModel.EVALUATION_VERSION`.
  - Not run: `test-fits`. No fit, and no quantity a fit's objective uses, changed; O3's fit uses its GVD, which is `dispersion_au` at order 2, unchanged.
- **Sources checked:** the three coefficients and their wavelength against their records (egan_2023 p. 17, egan_2026 pp. 4, 7).
- **Limitations:** the higher orders of `dispersion_au` for O3 remain gradients of the
  grid (the roadmap's item). The coefficients are not applied; §7.13's opt-in
  correction is deferred.
- **Open questions:** none new.
- **Handoff:** Phase 2c-i builds the equation-of-state object. OEP 4's Phase 1 can use
  `dpolarizability_domega` for the mixture's frequency derivative.

### Phase 2c, first part — The equation-of-state object

- **Scope:** the equation-of-state half of Phase 2c (§4.2, §6.5), 2026-09-30,
  implemented by Opus 5.5 (AI) for John Travers, on `feature/oep-0001-eos`, merged into
  `main`. Transport, the other half, follows as its own merge once its correlations'
  papers are checked (the author's decision of 2026-09-30).
- **Delivered:** `EquationOfState`, `equation_of_state`, `pressure`, with `phase_domain`,
  `saturation_pressure`, `speed_of_sound` and `compressibility`, in `optmat.gas.eos`;
  `number_density` and `compressibility` now delegate to it;
  `tests/test_equation_of_state.py`; the conventions page; `coolprop>=8` in the `gas`
  extra.
- **Deviations from the proposal:**
  - Phase 2c merges in two parts, so that OEP 4's Phase 2 is not held up by transport.
  - CoolProp is called directly, not through pyfluids: pyfluids does not expose the phase, the
    vapour pressure or the maximum pressure. The two agree bitwise (TC-758).
  - The domain ends 1e-5 below the vapour pressure (`SATURATION_MARGIN`), since the backend refuses to
    solve closer.
- **Decisions:** the round-trip tolerance, 1e-12 (`ROUND_TRIP_TOLERANCE`).
- **Register entries:** TC-758 (resolving TC-731).
- **Verification:** the full gate on the branch; §8's EOS row as TC-758 records.
- **Sources checked:** none new; the equations are CoolProp's, named by their keys.
  Their papers are curated with transport (the second part).
- **Limitations:** the equation's own uncertainty is not yet stated per gas; it comes with
  the second part's records.
- **Open questions:** none new.
- **Handoff:** OEP 4's Phase 2 builds its mixture equation of state on `EquationOfState`
  (its type) and `_abstract_state`. Phase 2d's air model uses CoolProp's pseudo-pure
  `Air`, which needs a species beside the `Gas` members.

### Phase 2d — Compositions, the air models and the air equation of state

- **Scope:** Phase 2d of §12 (§4.4, §6.6), 2026-09-30, implemented by Opus 5.5 (AI) for
  John Travers, on `feature/oep-0001-air`, merged into `main`.
- **Delivered:**
  - `composition.REDUCED`, and the two reduced names in `composition`;
  - `optmat.gas.air`: `AirModel`, `air_model`, `AIR_COMPOSITIONS` and `DEFAULT_AIR`;
  - `eos.air_equation_of_state`, and `EquationOfState.material` in place of a `Gas`;
  - the Raman response's named mixtures;
  - the source record `ciddor_1996`;
  - `tests/test_air.py`; the air models in the golden values; the models page.
- **Deviations from the proposal:**
  - `AirModel` is not a subclass of the gas `Model`, whose `gas` field is a member of
    `Gas`, which air is not. It carries the same interface, with its own identity.
  - §8 expects the reduced compositions to give solprop's recipes "to their rounding".
    O2 differs by 1.1e-4, beyond it, because solprop uses the older 20.95 % (TC-759). The test bounds the
    difference instead.
  - The Ciddor comparison uses his standard dry air at 15 °C and 1 atm, where no density
    equation of his enters, rather than his full moist-air procedure.
- **Decisions:** `"dry air"` names the default composition, for the air model as for the Raman
  response. `composition("dry air")` keeps its present meaning, the full table.
- **Register entries:** TC-759.
- **Verification:**
  - The full gate.
  - §8's air row: each model is the weighted sum of its constituents in value and derivative, to rounding.
  - The reduced compositions are CIPM-2007 renormalised, and against solprop within 1.2e-4.
  - The full composition agrees with Ciddor within 5e-8, and the omission and the density differences are as TC-759 records.
- **Sources checked:** Ciddor 1996's Eqs. (1)-(2) and Appendix A against the page image.
- **Limitations:**
  - The air model's pair terms, and so its index above a few bar, are OEP 4's.
  - The full composition is limited to 320–1400 nm by its molecules' present ranges, until OEP 3's mixture rule (its Q13) allows wider.
- **Open questions:** none new.
- **Handoff:** OEP 4's mixture optics generalises `AirModel` with pair terms. Its
  mixture equation of state replaces the pseudo-pure air for a labelled air fill.

### Phase 2c, second part — Transport of the dilute gases

- **Scope:** the transport half of Phase 2c (§4.3, §5, §6.5), 2026-09-30, implemented by
  Opus 5.5 (AI) for John Travers, on `feature/oep-0001-transport`, merged into `main`
  once each correlation's stated uncertainty was checked against its paper.
- **Delivered:**
  - `optmat.gas.transport`: `TransportModel`, `transport_model`, `viscosity`,
    `heat_capacity`, `molar_mass` and `mixture_viscosity` (Wilke), with
    `TRANSPORT_GASES` and `QUANTITIES`, for solprop's twelve gases and CO2;
  - `data/gas/transport.toml`, the curated choice of source per gas;
  - the quantities `dilute_viscosity` and `ideal_gas_heat_capacity` in a new family,
    `transport`, on {doc}`../../data`;
  - sixteen source records, seven with their printed dilute forms as formulas;
    `assael_1987`'s Table 2 as a dataset (`tools/extract/assael_1987.py`);
  - the WebBook's tables for Ne, Kr, Xe and N2O in `data/external/nist_webbook`;
  - `tests/test_transport.py`; the transport models in the golden values and the
    surface snapshot; the manual page {doc}`../../transport`.
- **Deviations from the proposal:**
  - Ne, Kr, Xe and N2O take Huber's dilute form, the models behind the WebBook's
    tables, not the tables at 0.1 MPa (§5, §6.5). It is the quantity §4.3 defines and
    covers each equation of state's temperatures; the tables stop at 200-400 K. The
    tables are kept as the check (TC-763).
  - §4.3's "a few parts in 10⁴" between one bar and the dilute limit understates it
    below room temperature: up to 1.7×10⁻³ at 200 K (N2), and 0.54 % for the WebBook's
    Xe. Every offset stays below its gas's stated uncertainty, so the rule stands
    (TC-762).
  - D2 (Q14): a primary was found, Assael et al. (1987). Its Table 2 is served,
    `provisional`, because its printed equations do not reproduce the table (TC-760).
  - SF6's correlation and equation of state are CoolProp's, unchecked, `provisional`:
    their papers are not in hand (TC-761). CH4 is `provisional` because no uncertainty is
    stated for its dilute term (TC-763).
- **Decisions:** $c_p^\circ$ keeps the CODATA convention of the number density: the
  equation's $c_p^\circ/R$ times CODATA's $R$, so the noble gases give $\tfrac52 R$
  exactly. Each viscosity's range is its paper's where the paper states one. CoolProp's
  viscosity is taken at 10⁻⁶ mol/m³ as its total, since CoolProp reports helium's as a
  single term (TC-763).
- **Register entries:** TC-760 to TC-763.
- **Verification:**
  - The full gate.
  - §8's transport row: each CoolProp correlation equals its paper's printed dilute
    form to 10⁻⁹. Huber's form reproduces his check values to 2.4×10⁻⁷. The noble
    gases' $c_p^\circ$ is $\tfrac52 R$ to 10⁻¹⁴.
  - Against solprop's sources (CoolProp at one bar, the WebBook at 0.1 MPa), optmat
    differs by the density term alone. D2's different source moves it by −0.41 % to
    +0.69 % over 200-600 K.
  - Wilke's rule reduces to the pure gas, does not depend on the order of the gases,
    and is 9×10⁻⁴ from Lemmon and Jacobsen's pseudo-pure air. Out-of-range calls are
    refused.
- **Sources checked:**
  - The dilute forms, constants, ranges and uncertainty statements of Arp et al. 1998,
    Lemmon and Jacobsen 2004, Muzny et al. 2013, Laesecke and Muzny 2017,
    Quiñones-Cisneros and Deiters 2006 and Huber 2018, in their text and tables.
  - Assael et al.'s Table 2 against the page image.
  - The ideal parts and molar masses of the ten equations of state against CoolProp:
    N2O's Eq. (10) directly, the others as each record says.
- **Limitations:**
  - Neon's equation of state (Thol et al.) is not obtainable, and CoolProp's helium
    equation is unpublished. Harvey et al.'s report stands for the latter.
  - Four papers are wanted: Assael et al. 1986 (H2), Quiñones-Cisneros et al. 2012 and
    Guder and Wagner 2009 (SF6), and Stewart et al. 1991 (O2's ideal part). Muzny et
    al.'s 2022 erratum is also wanted.
- **Open questions:** none new. Q14 is resolved, with D2 `provisional` (TC-760).
- **Handoff:** Phase 3's `solprop_transport.json` records solprop's values against these
  models. OEP 4's mixture work may take the molar masses and heat capacities from here.

### Phase 2a — Explicit parameters and overrides

- **Scope:** Phase 2a of §12 (§6.3, §7.4, §7.5), 2026-09-30, implemented by Opus 5.5
  (AI) for John Travers, on `feature/oep-0001-overrides`, merged into `main`.
- **Delivered:**
  - `Identified.with_parameters`, with `INPUTS` and `DERIVED` per class, and
    `modified` and `base` on every model and in its identity; `identity.describe`;
  - `from_parameters` on both rate models and on `OscillatorModel` (SI);
    `GasRamanResponse.from_lines`; `FormulaModel.sellmeier` and
    `TabulatedModel.from_arrays` for the solids;
  - the inputs of every class of the table on {doc}`../../identity`;
  - custom species as plain names in the rate and gas dispersion models;
  - `tests/test_overrides.py`; the identity page's section; the surface snapshot.
- **Deviations from the proposal:**
  - A dispersion model's overrides name its fields in atomic units (§7.5), while
    `OscillatorModel.from_parameters` takes SI (strengths, pole energies in J, an
    offset in C m²/V; OEP 2's form).
  - A fitted model's term count is fixed except through `infrared`.
  - `MolecularRateModel` keeps no structure coefficients, so an override that does not
    name them passes the curated orientation average through (`structure_factors`).
  - The equation-of-state, transport, air, hyperpolarizability, Kerr and
    measured-absorption models declare no inputs: §6.3 gives no rules for them.
  - The Raman density range needs each band's strongest line's Doppler width. The
    response stores only the Dicke caps, so `from_lines` computes it from the
    species' mass (curated, or `mass=`).
  - `with_ionization_energy` is an override, so its result is `user`.
- **Decisions:** a custom species takes the ideal-gas density and no virial ratio;
  `modified` and `base` are outside the digest, as the status is (TC-764).
- **Register entries:** TC-764.
- **Verification:**
  - The full gate.
  - §8's override row:
    - `from_parameters` with the curated inputs reproduces every curated rate and
      molecular rate, digest for digest, and `from_lines` every tested Raman response.
    - The curated builders go through them, and every golden value and digest is
      unchanged.
    - For each input tested, `with_parameters(x=v)` equals the constructor of the
      changed inputs in digest and numbers, derived values included. That covers the
      structure factors, the Hartree coefficient, and the Raman groups and range
      after a width law or scale.
    - Derived fields are refused, a Dicke term without a cap is refused, and
      `modified`, `base` and the status are set.
- **Sources checked:** none new.
- **Limitations:** a custom species cannot use the curated equations of state or
  virial data; retyping `Gas` across the package is TC-317's follow-up.
- **Open questions:** none new.
- **Handoff:** OEP 2's compact models can be built with `OscillatorModel.from_parameters`.
  Phases 2e (Raman groups) and 2f (rate energies) build on the rebuild paths here.

### Phase 2c, follow-up — The transport papers

- **Scope:** the papers the second part asked for, curated after its merge, 2026-09-30,
  by Opus 5.5 (AI) for John Travers, on `feature/oep-0001-transport-followup`, merged
  into `main`.
- **Delivered:** the records `assael_1986`, `quinones_cisneros_2012`, `guder_2009` and
  `stewart_1991`; `TransportModel.heat_capacity_range` and the choices file's
  `heat_capacity_range_k`; tests reproducing each paper's printed form.
- **Deviations from the proposal:** none. Q14 is now fully settled: D2's primary is
  Assael et al. (1987), and its table is explained by their 1986 hydrogen correlation.
- **Decisions:** D2 and SF6 are `validated`; O2's heat capacity stops at 1700 K, where
  CoolProp's ideal part leaves the stated accuracy of its source (TC-765).
- **Register entries:** TC-765 (resolving TC-760 and TC-761).
- **Verification:** the full gate; the golden values re-recorded with no value
  changed. Only the transport digests moved, since the model gained a field.
- **Sources checked:**
  - Assael et al. 1986, Eqs. (1), (3), (13) and Table 2.
  - Quiñones-Cisneros et al. 2012, Eq. (6), Table 2 and the uncertainty statement.
  - Guder and Wagner 2009, Eq. (15) and Table 21 (against the image).
  - Stewart et al. 1991, Eqs. (11) and (12).
  - Muzny et al.'s 2022 erratum.
- **Limitations:** Guder and Wagner's uncertified OCR tables are not used; neon's
  equation-of-state paper remains unobtainable.
- **Open questions:** none new.
- **Handoff:** none.

### Phase 2e — The Raman groups

- **Scope:** Phase 2e of §12 (§4.5, §6.7), 2026-09-30, implemented by Opus 5.5 (AI) for
  John Travers, on `feature/oep-0001-raman-groups`, merged into `main`.
- **Delivered:** `RamanGroup` and `GasRamanResponse.groups()`, `line_id` and
  `check_densities`; short band keys in `raman.toml`; the line-identifier snapshot and
  its tool; `tests/test_raman_groups.py`; the manual page's sections.
- **Deviations from the proposal:** the identifiers use a short curated band key (`v1`,
  `nu3`, `fermi_1388`), since the curated band names are long and one contains "/";
  the response's `band` is that key, so the Raman digests moved (TC-766).
- **Decisions:** the group view never warns; `check_densities` warns for declared
  densities outside `density_range`.
- **Register entries:** TC-766.
- **Verification:**
  - The full gate.
  - §8's group row:
    - The members' width laws and shifts are equal to 0 ulp for every gas.
    - The groups sum to `response` to roundoff, over densities spanning the range.
    - The identifiers are unique and match the snapshot.
    - `decay_rate(0)` is the limit, and the rate tends to it.
    - The copolarised and cross-polarised elements, and solprop's weights, are
      reproduced line by line.
    - Only H2's Q group shifts.
- **Sources checked:** none new.
- **Limitations:** `response` and `spectrum` are not vectorised over densities
  (optional in §6.7).
- **Open questions:** none new.
- **Handoff:** solprop's Raman leaf can take the groups directly (SEP 3).

### Phase 2f — The rate energies

- **Scope:** Phase 2f of §12 (§4.6, §6.8), 2026-09-30, implemented by Opus 5.5 (AI) for
  John Travers, on `feature/oep-0001-rate-energies`, merged into `main`.
- **Delivered:** `channel_rates`, `energy_loss_rate`, `energy_per_event` and
  `energy_kind` on both rate models, with the linear-polarisation statement in their
  docstrings; `tests/test_rate_energies.py`; the rates page's section.
- **Deviations from the proposal:** none.
- **Decisions:** an atom's `energy_kind` is its curated kind (`adiabatic`), or `not
  stated` for a custom species or an overridden energy (TC-767).
- **Register entries:** TC-767.
- **Verification:**
  - The full gate.
  - §8's rate-energy row: the channel rates sum to the rate (1e-12), and the loss
    rate is the channel sum of energy times rate, 0 at zero field.
  - The energy per event is their ratio, the lowest channel's energy at zero field
    and where the rates underflow, and a single channel's own energy.
  - Argon's 15.797 and 15.814 eV of §4.6 are reproduced.
- **Sources checked:** none new.
- **Limitations:** none new.
- **Open questions:** none new.
- **Handoff:** solprop's ionisation leaf can integrate `energy_loss_rate` beside the
  rate.

### Phase 2g — The silica bundle

- **Scope:** Phase 2g of §12 (§4.7, §6.9, §7.10), 2026-10-01, implemented by Opus 5.5 (AI)
  for John Travers, on `feature/oep-0001-silica-bundle`, merged into `main`.
- **Delivered:** `optmat.solid.bundle` with `SilicaNonlinearBundle`, `silica_bundle`,
  `silica_bundles` and `GainCheck`; the curated tensor default; the n0 = 1 rule in
  `raman_coefficient`'s docstring; `tests/test_silica_bundle.py`; the silica page's
  section; the bundles in the golden values and the surface snapshot.
- **Deviations from the proposal:** the bundle's module is `optmat.solid.bundle`, a
  new module beside the Kerr fit and the responses (it imports both). It adds
  `n2_raman` beside §6.9's `n2` and its electronic part.
- **Decisions:** the per-response bundles of the default and tensor responses are
  `validated` like those two; Blow and Wood's is `comparison` (TC-768).
- **Register entries:** TC-768.
- **Verification:**
  - The full gate.
  - §8's silica row: each bundle's total is the Kerr fit's n2, and its parts add up
    to it; `paired` means the fit's f_R; every bundle's gain check lies within Stolen
    and Ippen's 10 %.
  - The gains of §4.7's table are reproduced: −9.7 % and −7.9 %.
- **Sources checked:** none new.
- **Limitations:** the tensor pairing's option (c) of §7.10, a second Kerr fit with Lin
  and Agrawal's shape, is not done.
- **Open questions:** none new.
- **Handoff:** solprop's silica leaf takes `silica_bundle("default")` for a scalar run
  and `silica_bundle("tensor")` for a vector one.

### Phase 2h — The solids

- **Scope:** Phase 2h of §12 (§6.10), 2026-10-01, implemented by Opus 5.5 (AI) for John
  Travers, on `feature/oep-0001-solid-tables`, merged into `main`.
- **Delivered:** `TabulatedModel.out_of_range` (`"hold"` or `"raise"`) and
  `with_out_of_range`; the capability notes on the metals and crystals pages;
  `tests/test_solid_tables.py`. Identity and status on every solid model came with
  Phase 1.
- **Deviations from the proposal:** the rule is set by `with_out_of_range`, not by an
  override: it changes no returned number, so the model keeps its status (TC-769).
- **Decisions:** a NaN wavelength is not refused.
- **Register entries:** TC-769.
- **Verification:** the full gate; the golden values re-recorded with no value changed.
- **Sources checked:** none new.
- **Limitations:** none new.
- **Open questions:** none new.
- **Handoff:** Phase 2 is complete. Phases 3-5 follow OEP 3's refits (§12).

### Phase 3 — Parity fixtures

- **Scope:** Phase 3 of §12 (§6.11), 2026-10-01, implemented by Opus 5.5 (AI) for John
  Travers, on `feature/oep-0001-parity`, merged into `main`.
- **Delivered:**
  - every `tests/data/solprop_*.json` re-recorded at solprop main d57741490dbc, the
    commit inside each;
  - the new `solprop_transport.json` (viscosity, heat-capacity ratio and molar mass over
    each gas's range) and `solprop_eos.json` (number density and the pressure round
    trip), with their tools;
  - `tests/test_solprop_baseline.py`;
  - the comparison table on {doc}`../../transport`.
- **Deviations from the proposal:**
  - solprop exposes the heat-capacity ratio rather than $c_p$, so the fixture records
    $\gamma$, compared with optmat's $c_p^\circ/(c_p^\circ - R)$.
  - Air is in the density fixture (solprop's pseudo-pure air) but not in the transport
    one (solprop has no transport entry for it).
- **Decisions made:** none beyond the fixtures' grids (nine temperatures across each
  range; 0.1-30 bar at 273.15 and 294 K).
- **Register entries:** TC-787.
- **Verification:**
  - The full gate.
  - The seven re-recorded fixtures' values are identical to the previous recordings.
  - The densities agree to the gas-constant convention within 3e-10.
  - The transport table is recomputed.
- **Sources checked:** solprop's transport and equation-of-state modules (their
  sources, ranges and conventions) against its code at d577414.
- **Findings:**
  - solprop returns a liquid density for SF6 at 30 bar.
  - The transport page had the sign of the D2 difference reversed; it is now
    corrected.
- **Limitations:** the fixtures are solprop's values at one commit. They are re-recorded
  when solprop or the OEP 3 refits change them.
- **Open questions:** none new.
- **Handoff:**
  - Phase 4 writes the page for propagation codes.
  - When OEP 3's refits land, compare the fixtures holding CO, CO2, CH4 and N2O again
    (§6.11).

### Phase 4 — The page for propagation codes

- **Scope:** Phase 4 of §12 (§9), 2026-10-01, implemented by Opus 5.5 (AI) for John
  Travers, on `feature/oep-0001-propagation-page`, merged into `main`.
- **Delivered:** {doc}`../../propagation`, with every item of §9's list. Its snippets run in
  the docs test, and it links the pages that hold the detail ({doc}`../../identity`,
  {doc}`../../transport`, {doc}`../../mixtures`, {doc}`../../raman_response`,
  {doc}`../../ionisation_rates`, {doc}`../../silica`).
- **Deviations from the proposal:** the snippets run in `tests/test_docs_snippets.py`,
  as every page's do, rather than as doctests.
- **Decisions made:** none.
- **Register entries:** TC-788.
- **Verification:** the full gate; the page's snippets.
- **Sources checked:** none; the page documents the interface.
- **Limitations:** the page states the interface as it is now. OEP 3's refits will
  change the molecules' intervals, and the page names where to read them.
- **Open questions:** none new.
- **Handoff:** Phase 5 adds the clean-install test.

### Phase 5 — Release

- **Scope:** Phase 5 of §12 (§6.12), 2026-10-01, implemented by Opus 5.5 (AI) for John
  Travers, on `feature/oep-0001-release`, merged into `main`. The tag is left to the
  author.
- **Delivered:**
  - `tests/test_clean_install.py`, marked `release`, run by `poe test-release` and left
    out of the default lane.
  - It builds the wheel and the sdist, installs each with the `gas` extra offline into a
    fresh environment away from the source tree and without the paper library, and
    runs the quickstart's and the propagation page's snippets and the whole public
    surface there.
  - It also checks the artifacts' contents.
- **Deviations from the proposal:** the clean environment installs from uv's cache
  rather than with no package source at all. That is offline, but needs the
  dependencies cached.
- **Decisions made:** the sdist excludes `/.claude` and `/.codex`.
- **Register entries:** TC-789.
- **Verification:**
  - The release lane passes for both artifacts.
  - The full gate passes.
- **Findings:** the sdist held 265 MB of agents' worktrees before the exclusion
  (TC-789).
- **Limitations:**
  - The test runs on the developer's platform only; the wheel is pure Python.
  - The golden digests' platform dependence is Q32's.
- **Open questions:** Q31 and Q32 remain from Phase 1. They are moved to the roadmap's
  *Infrastructure* (2026-10-01).
- **Every phase is merged and reported.** The Final status is the author's to approve
  (OEP 0), with the tag. Before it:
  - the open questions are to be settled or moved to the roadmap;
  - each page OEP 1 produced takes its closing link to the OEP;
  - when OEP 3's refits land, the parity fixtures of the refitted gases are compared
    again (§6.11).

## Part III — History

- 2026-09-29: Drafted. Its recommendations were discussed before the OEP process
  existed. Here they are proposals, to be settled in review, and no earlier approval is
  carried over.
- 2026-09-29: GPT-6 (AI) reviewed all of Part I against the current optmat evaluators,
  manual and register, and solprop SEP 3. Recorded unresolved findings as Q18--Q27;
  the recommendation is to revise before acceptance. No implementation or scientific
  default was changed. The review checked existing code and curated records, with
  numerical probes of overrides, low-density widths and the silica gain; it did not
  freshly verify the underlying papers or close Q17.
- 2026-09-30: Q28 added: the polarizability's own temperature dependence, after the
  molecular fits took a reference temperature (TC-730). No recommendation of §7 changed.
- 2026-09-30: Opus 5.5 (AI) reviewed all of Part I against the code at `ea0dcaf`, the
  register, solprop main at `8a4093f` and SEP 3, and resolved Q13, Q15, Q16 and Q18–Q27
  in the text. It ran the checks that could be run: the silica gain of every bundle
  (the paired Lin–Agrawal variant is 17.8 % below the measured gain), the equation of
  state across phases (the forward law returns liquid densities past saturation), the
  low-density limit of every curated width law, the rate channels of the noble gases,
  O3's measured-absorption derivative, and the air density rules. §7.3, §7.10 and the
  new §7.13 changed their recommendations; Q29 was added. Left for the maintainers: the
  confirmations of §7, Q10, Q17's check of Wilke's paper, and Q28; D2's viscosity is
  deferred to Phase 2c. The underlying papers were not freshly verified, except Huber's
  report for §5 and the DOI records of the † references.
- 2026-09-30: J. C. Travers accepted the review's recommendations: the §7
  recommendations confirmed; Q10 decided for the unpaired Lin–Agrawal tensor default
  (TC-734); Q28 decided for metadata now and the correction later; Q29 decided in
  SEP 3 for the pseudo-pure air equation. The review's findings entered the register as
  TC-731 to TC-734. Open before acceptance: Q17's check of Wilke's paper, and the
  human review of Part I.
- 2026-09-30: Q30 added, with the surface OEP 3 needs (qualified intervals, their
  identity, the infrared terms' temperature, the validity of air, the infrared terms as
  one override input, `OscillatorModel.from_parameters`, and the fixtures after OEP 3's
  refits). No recommendation of §7 changed, and no number or warning changes.
- 2026-09-30: John Travers (human) reviewed Part I. With the AI reviews, the review
  status is `Reviewed`. Still open before acceptance: Q17's check of Wilke's paper.
- 2026-09-30: Aligned with OEP 3's revision for its review. OEP 3's representation
  uncertainty is now a representation bound (§4.8). An intersection of intervals widens
  only under OEP 3's three checks (§4.4, Q30). The gaps between OEP 3's intervals carry
  a reason (§4.1). No recommendation of §7 changed, and no number or warning changes.
- 2026-09-30: Q17 resolved: Wilke's paper and Lemmon et al.'s air equation curated
  (TC-737, TC-741) and checked against the pages. §5 names the primary paper of
  CoolProp's pseudo-pure air.
- 2026-09-30: Accepted by John Travers. Part I is frozen.
- 2026-09-30: Implementation started. Phase 1 (surface and identity) merged, with its
  report in Part II (TC-745, TC-746).
- 2026-09-30: Phase 2b (the polarizability derivative and the reference temperature)
  merged, with its report in Part II (TC-757).
- 2026-09-30: Phase 2c's first part (the equation-of-state object) merged, with its
  report in Part II (TC-758).
- 2026-09-30: Phase 2d (compositions, the air models and the air equation of state)
  merged, with its report in Part II (TC-759).
- 2026-09-30: Phase 2c's second part (transport of the dilute gases) merged, with its
  report in Part II (TC-760 to TC-763).
- 2026-09-30: Phase 2a (explicit parameters and overrides) merged, with its report in
  Part II (TC-764).
- 2026-09-30: Phase 2c's follow-up (the transport papers) merged, with its report in
  Part II (TC-765).
- 2026-09-30: Phase 2e (the Raman groups) merged, with its report in Part II (TC-766).
- 2026-09-30: Phase 2f (the rate energies) merged, with its report in Part II (TC-767).
- 2026-10-01: Phase 2g (the silica bundle) merged, with its report in Part II (TC-768).
- 2026-10-01: Phase 2h (the solids) merged, with its report in Part II (TC-769). Phase 2
  is complete.
- 2026-10-01: Every phase merged and reported; Q31 and Q32 moved to the roadmap. Status
  Final, approved by the author, and the OEP archived.
