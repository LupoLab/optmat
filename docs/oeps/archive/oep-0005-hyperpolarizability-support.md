# OEP 5 — A public record of where the hyperpolarizability fits are supported

| | |
|---|---|
| **OEP** | 5 |
| **Title** | A public record of where the hyperpolarizability fits are supported |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5; GPT-6 (review revisions, 2026-10-07) |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | John Travers (human); GPT-6 (AI) |
| **Type** | Feature |
| **Created** | 2026-10-07 |
| **Resolution** | 2026-10-07: accepted by the author, John Travers, after the human and AI reviews; merged directly to `main` |
| **Requires** | [OEP 1](oep-0001-propagation-interface.md): model identity, digests and the declared public surface |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phase 1 merged from `feature/oep-0005-support-record` |
| **Manual pages** | {doc}`../../nonlinear` (*The support record*); {doc}`../../identity` (*Records about a model*); {doc}`../../propagation` (*Where a model holds*) |

## Part I — Proposal

### 1. Abstract

optmat's fitted second hyperpolarizability of each gas,
$\gamma = \gamma_0[1 + A W_2 + B W_2^2 + B' W_4]$, is about to be used by solprop as a
full frequency-dependent electronic $\chi^{(3)}$: every four-wave quartet
$(-\nu_\sigma; \nu_1, \nu_2, \nu_3)$ is evaluated through filters that act on each
field, that is on each leg separately (solprop's SEP 6). Broadband propagation then
meets quartets outside the fitted data, and solprop needs, per gas, where the fit is
supported and where it cannot be trusted, without copying numbers from optmat's manual.
This OEP proposes a documented, tested, read-only record of material facts on
`HyperpolarizabilityModel`: the largest fitted $W_2$, single leg and pair sum with the
values that set them, and the gas's resonance limits with the curated datasets they come
from, plus two quantities derived from them by definition. What a propagation code does
with them is its own modelling decision (TC-868); the record carries no rule. It does
not change any evaluated gamma or any stored fit. It succeeds if solprop's SEP 6 can
build its domain policy from the record alone, and if every number in it is traceable
and tested.

### 2. Motivation and scope

**The need.** solprop's SEP 6 (*Frequency-dependent electronic third-order response*,
Draft) asks, in its open question Q2, what public optmat record expresses a qualified
quartet domain, its resonance exclusions and its uncertainty. Its demonstration cases are
deep-ultraviolet four-wave mixing in argon (400 + 400 - 800 nm → 266.7 nm) and
vacuum-ultraviolet resonant dispersive waves in helium (100-150 nm from an 800 nm pump).
solprop's fallback today is a constant Kerr coefficient at the carrier. That is also an
extrapolation. For the stored polynomials, increasing every leg magnitude above the
carrier's increases gamma, but a mixed quartet need not lie above the carrier value,
and the physical response has no such general ordering. Comparing a proposed
extrapolation with the carrier is solprop's task; this record does not certify either.

**What optmat publishes today is not enough.** A model carries its declared range
(`range_nm`, a fundamental wavelength) and the largest fitted $W_2$
(`record["max_w2_fitted_cm-2"]`), and warns beyond 1.1 times the latter. Neither says
which quartets the fit supports:

- capping every leg at the Kerr diagonal of the largest fitted $W_2$ (374 nm for argon)
  changes the argon quartet above, which lies within the fitted $W_2$ (0.98 of it), by
  -9.3 %, against -4.3 % for a 400 nm carrier;
- the single legs of the data reach 1.6 to 1.7 times further out in wavenumber than that
  diagonal, because the visible ESHG and the calculated third harmonics put their
  harmonic legs at 229-244 nm;
- the resonances that end any extrapolation are a property of the gas, not of the fit,
  and are not recorded with it.

**In scope:**
- a read-only record of each gas's support and resonance limits, with two derived
  quantities, on the model, with its own digest;
- the curated data it needs (the resonance levels) and the stored, derived numbers;
- tests, docstrings, a manual section and the register entries.

**Out of scope:**
- any change to `gamma_au` or the stored fits, and any refit;
- a resonant model of any gas (solprop's resonant helium response is SEP 12's);
- the curation of the polyatomic molecules' two-photon states and absorption onsets
  (data work, requested on the roadmap; the record says where they are missing);
- any rule for using the facts: optmat publishes facts, and what to admit, cap or label
  is solprop's decision (SEP 6, section 4.10; TC-868).

### 3. Context

- **Manual:** {doc}`../../nonlinear` (the model, its fits, the time-domain form and its
  warning beyond the fitted $W_2$) and {doc}`../../nonlinear_accuracy` (how well each fit is
  known). The analysis this proposal rests on is there: the legs of the fitted data,
  the resonance limits of every gas, the monotonicity of the fits and the quantities
  derived from them ({doc}`../../nonlinear`, *Where the fits are supported*), and the
  evaluation of
  per-leg bounds against the data, a Kerr carrier and solprop's demonstration quartets
  ({doc}`../../nonlinear_accuracy`, *Beyond the fitted data*). The development tool
  `tools/fit/nonlinear_support.py` computes those pages' tables, and
  `tests/test_nonlinear_support_page.py` checks them.
- **Modules:** `optmat.gas.nonlinear` (`HyperpolarizabilityModel`, `frequency_sums`,
  `W2_MARGIN`, `hyperpolarizability_models`); `optmat.identity` (`Identified`).
- **Register:** TC-313 to TC-316 (the model, its fits and its warning beyond the fitted
  $W_2$; the 650 nm window on calculated third harmonics); TC-745 (identity of the
  models); TC-443 and TC-444 (NIST ASD energies and resonance lines); TC-854 (xenon's
  negative n2 above its two-photon resonance); TC-862 and TC-863 (the resonance levels
  of the atoms and the electronic states of the diatomic molecules); TC-864 (the legs of
  the fitted data, and the correction of the calculated third harmonics' reach to
  231.4 nm); TC-865 (the resonance limits, open for CO2, SF6 and every polyatomic
  two-photon state); TC-866 (monotonicity); TC-867 (the evaluation of per-leg bounds,
  kept as evidence); TC-868 (optmat states facts; their use is the propagation
  code's).
- **Roadmap:** the *solprop* item, "the process-aware Kerr and third-harmonic
  coefficients, through the model's exact time-domain form".
- **Prior OEPs:** OEP 1 (the propagation-code interface: identity, digests, overrides
  and the declared surface).

### 4. Theory

Wavenumbers here are signed inverse vacuum wavelengths in m⁻¹, not angular
wavenumbers: $\omega_k = 2\pi c s_k$. Write $W_2 = \sum_k s_k^2$ in m⁻² and
$W_4 = \sum_k s_k^4$ in m⁻⁴. Thus this section's $W_2$ is SEP 6's $S_2$, or $10^4$
times the numerical value returned by `frequency_sums` with cm⁻¹ inputs; $W_4$
converts by $10^8$. In the polynomial, the corresponding SI coefficients are the
stored numerical $A$, $B$, $B'$ times $10^{-14}$ m², $10^{-28}$ m⁴ and $10^{-28}$ m⁴.
The symbol $E_i$ below means a transition energy divided by $hc$, in m⁻¹.

**Legs.** Write each quartet as four signed wavenumbers $s_k$ that sum to zero, the output
with a minus sign: the Kerr effect is $(-\nu, \nu, \nu, -\nu)$, the third harmonic
$(-3\nu, \nu, \nu, \nu)$, ESHG $(-2\nu, \nu, \nu, 0)$, the dc Kerr effect
$(-\nu, \nu, 0, 0)$, CARS $(-(2\nu - \nu_s), \nu, \nu, -\nu_s)$ and the pump-probe index
$(-\nu_r, \nu_r, \nu, -\nu)$. Then $W_2 = \sum_k s_k^2$ and $W_4 = \sum_k s_k^4$, as
`frequency_sums` gives them. Since each leg is minus the sum of the other three, the
Cauchy-Schwarz inequality gives $s_j^2 \le 3(W_2 - s_j^2)$, so

$$
\max_j |s_j| \le \sqrt{3 W_2 / 4},
$$

with equality for the third-harmonic output. No quartet within the largest fitted $W_2$
has a leg beyond $\nu^\ast = \sqrt{3 W_{2,\max}/4}$.

**Resonances.** This is the electronic response in the electric-dipole approximation,
with the curated ground-state and spin selection assumptions (TC-862, TC-865).
Rotational and vibrational Raman resonances, linewidths, thermally populated states
and weak forbidden transitions are not qualified by these limits. Being below a level
does not establish a specified accuracy or negligible absorption.

In the connected sum over excited states for gamma (Shelton & Rice [1], eq. 14 on
p. 8, and eq. 17 on p. 10 for the static limit), terms contain four
dipole matrix elements
$\langle g|d|m\rangle\langle m|d|n\rangle\langle n|d|p\rangle\langle p|d|g\rangle$ and
three energy denominators. For one input ordering they are
$\Omega_{mg} - (\omega_1+\omega_2+\omega_3)$,
$\Omega_{ng} - (\omega_1+\omega_2)$ and $\Omega_{pg} - \omega_1$, where
$\Omega_{ag}=(\mathcal E_a-\mathcal E_g)/\hbar$ is a positive excitation frequency.
Permutations and counter-rotating terms give the other signs. For polar molecules
the internal dipoles are ground-dipole-subtracted; the full expression also has
disconnected subtraction terms. This is a description of possible poles, not a
complete expression for gamma: transition amplitudes can vanish. The output leg
accounts for three-photon resonances of the inputs. Within the stated approximation:

- $E_1$, the lowest included electronic level dipole-coupled to $g$, gives the
  below-level condition $|s_j| < E_1$;
- $E_2$, the lowest included excited electronic level two dipole steps from $g$, gives
  the below-level condition for every pair sum:
  $|s_a + s_b| < E_2$, $a\ne b$ (two-photon absorption for legs of one sign, an
  electronic Raman resonance for opposite signs);
- the ionisation energy, reached by a leg or a pair, opens a continuum and makes gamma
  complex. It lies above the selected discrete levels where both are curated;
  missing levels do not establish that ordering for a molecule.

For an atom with a $J = 0$ ground level, $E_1$ is the lowest odd $J = 1$ level and $E_2$
the lowest included excited even level with $J \le 2$. The two-photon amplitude to
an even level is a tensor in the two photons' polarisations $\mathbf e_a, \mathbf e_b$;
$J = 0$ and 2 are
reached through its symmetric part, $J = 1$ through its antisymmetric part, which is
proportional to $\mathbf e_a \times \mathbf e_b$ and to
$1/(\omega_{mg} - \omega_a) - 1/(\omega_{mg} - \omega_b)$, and so vanishes for one linear
polarisation (a scalar field) and for equal frequencies. The record therefore carries
two pair limits, for a scalar field and for any polarisation. The equal-frequency
cancellation is explicit in DeMille et al. [9], eqs. 1-4. Molecules follow the same
logic with their own symmetries; without inversion symmetry, one- and two-photon
transitions need not be mutually exclusive. Each proposed assignment still requires
its own selection-rule justification. Two polarisation-dependent level lists do not
qualify a scalar fitted gamma as an arbitrary tensor response.

**A sufficient condition on all four legs.** If every leg, including the generated
output, has magnitude at most $K$, then $|s_a+s_b|\le 2K$. Thus
$K<\min(E_1,E_2/2)$ keeps all leg and pair magnitudes below the selected electronic
levels. A bound on the input band alone does not bound newly generated outputs;
the consumer must include those outputs in the set being assessed.

For a scalar field, the SPM quartet $(-\nu,\nu,\nu,-\nu)$ has pair magnitude
$2\nu$, and can attain the scalar pair threshold when the transition is allowed.
For general polarisation, $\min(E_1,E_{2,\mathrm{any}}/2)$ remains a sufficient
below-level bound, but need not be an SPM resonance: an atomic even $J=1$ transition
vanishes at equal photon frequencies [9]. Exceeding either bound does not prove that
an occupied quartet is resonant, particularly for a sparse spectrum. The record gives
neither linewidths nor tensor amplitudes. A harmonic leg near $E_2/2$ can pair only
with smaller legs, so this band bound is not a necessary condition on each quartet.

**Monotonicity.** For one leg, $\partial(\gamma/\gamma_0)/\partial(s_j^2) =
A + 2BW_2 + 2B's_j^2$. With $\gamma_0>0$, $A > 0$, $B \ge 0$, $B + B' > 0$
(true of every stored fit)
and $s_j^2 \le W_2$, it is at least $A$ at any $W_2$. This is an algebraic property
of the polynomial on nonnegative squared legs, including capped tuples that no longer
sum to zero; it is not a physical continuation of a frequency-conserving quartet.
The physical gamma has no such theorem: the static sum over
states has a negative term, and a negative static gamma is possible in principle;
Shelton & Rice reported in 1994 that no such physical case was then known [1].
It rises toward the first resonance in the curated calculations, and above a
two-photon resonance it can change sign (xenon at 248 nm [2]).

**Facts and their use.** The record states facts ($W_{2,\max}$, the largest fitted leg
and pair, $E_1$, $E_2$, the ionisation energy) and, labelled as derived, $\nu^\ast$ and
$\min(E_1, E_2/2)$. How to use them is the propagation code's decision (TC-868). The
manual's accuracy page gives the evidence such a decision can draw on (TC-867): a leg
filtered to $\min(|s_j|, \nu^\ast)$ leaves every fitted value and every quartet within
the fitted $W_2$ unchanged; beyond the data the filtered value is at or below the
polynomial.
The polynomial generally underestimates the tabulated calculations farther beyond
the fitted windows, but not every point: the argon THG row of *Beyond the fitted data*
starts at a model/value ratio of 1.001. These comparisons cover particular processes
and gases and include fitted scale factors; they establish no universal lower bound
on the physical response of arbitrary mixed quartets or an error bar on extrapolation.
solprop's rule (SEP 6, section 4.10) caps every leg at $\nu^\ast$ and labels the band
supported, bounded or best effort.

### 5. Primary sources and data

- **NIST Atomic Spectra Database** (Kramida, Ralchenko, Reader and NIST ASD Team,
  doi:10.18434/T4W30F), curated as `nist_asd_2026` (TC-443): the ionisation energies and
  resonance lines, and the neutral levels of He to Xe stored as returned (TC-862). In
  hand.
- **Huber & Herzberg's constants of diatomic molecules** [3], through the NIST Chemistry
  WebBook (doi:10.18434/T4D303), curated as `nist_webbook_2026`: the electronic states of
  H2, D2, N2, O2 and CO, stored as returned (TC-863). In hand; the book itself is not
  curated (marked).
- **The polyatomic molecules' absorption:** Chan, Cooper & Brion's N2O (`chan_1994`);
  Ityaksov, Linnartz & Ubachs's CH4 and SF6, no detected absorption at their sampled
  wavelengths down to 197.7 nm
  (`ityaksov_2008`); Serdyuchenko et al.'s O3 (`serdyuchenko_2014`). In hand. Not curated
  and requested: CO2's and SF6's absorption onsets, and every polyatomic two-photon
  state.
- **Ionisation energies:** the datasets of `optmat.gas.ionisation`, and Mahla et al.'s
  measured thresholds for O3 (`mahla_2025`).
- **Vocabulary:** `level_energy` (atoms and ions) and `electronic_state_energy`
  (a molecule's excited state, TC-863); both exist.
- **Validators:** the calculated values outside the fitted windows (Bishop & Pipin's
  helium to 227.8 nm, Hättig & Jørgensen's Ne, Ar and N2 third harmonics, Bishop et al.'s
  H2 and D2) and the validators beyond the fitted $W_2$ test extrapolation on those
  processes. They do not qualify arbitrary quartets.
- **Uncertainty:** preserve each source's stated uncertainty and confidence level,
  or explicitly say it is unstated. TC-862 records missing uncertainties for Ar, Kr
  and Xe; printed precision is not an uncertainty. Molecular $T_e$ (potential minima)
  and $T_0$ (vibrational origins) describe different quantities: their difference is
  not an error bar. Record which was used, the chain of band origins and any
  approximation. Propagation of uncertainty through a chain must state its treatment
  of correlations; absent source information must not become zero uncertainty.
- **Absorption evidence:** CH4 and SF6 have no detected absorption at the sampled
  wavelengths of [6], not a proven absence of all lower electronic resonances.
  N2O's 5.70 eV is an integration boundary of [7]. O3's 1046.8 nm is the longest
  wavelength of a curated table, not an absorption onset; [8] reports measurements
  over 213-1100 nm. Retain these as explicitly qualified evidence, not precise levels.

### 6. Proposed interface

A frozen dataclass `HyperpolarizabilitySupport` in `optmat.gas.nonlinear`, returned by a
new read-only property `HyperpolarizabilityModel.support`, and a small
`SupportLimit` for each limit. Wavenumbers in m⁻¹; wavelengths derivable. Nested
records and sequences are immutable too. The property is optional for a model
without a matching curated support record (see *Identity* below).

| Name | Meaning | Units |
|---|---|---|
| `max_w2` | the largest fitted $W_2$ (today's stored value) | m⁻² |
| `max_w2_witnesses` | the fitted points attaining `max_w2`, with row-level provenance | tuple |
| `leg_max` | the largest fitted single leg, a `SupportLimit` naming its dataset, process and fundamental | m⁻¹ |
| `pair_max` | the largest fitted pair sum, likewise | m⁻¹ |
| `e1` | the lowest included dipole-coupled level, or qualified absorption evidence (`SupportLimit`, or None) | m⁻¹ |
| `e2`, `e2_any` | the lowest included two-photon level for a scalar field and for any polarisation (or None) | m⁻¹ |
| `ionisation` | the ionisation energy, with its dataset | m⁻¹ |
| `leg_max_in_w2` | derived: $\sqrt{3 W_{2,\max}/4}$, the largest leg of any quartet within the fitted $W_2$ | m⁻¹ |
| `band_edge(polarisation="linear")` | derived: $\min(E_1, E_2/2)$, the sufficient bound of section 4; None if either selected limit is missing or is only absorption/bound evidence | m⁻¹ |
| `model_digest` | the evaluated model whose fitted support is described | str |
| `digest` | the record's own identity (OEP 1) | str |

`digest` is a read-only convenience property for `identity.digest`, using OEP 1's
canonical identity machinery; it is not a second hashing convention.

Each `SupportLimit` carries `wavenumber`, a uniform `datasets` tuple, the selected
rows/columns or source-text location, the derivation and its register citations.
Its machine-readable `kind` distinguishes `fitted_extent`, `level`, `band_origin`,
`lower_bound`, `upper_bound` and `absorption_evidence`. Bound direction is always in
wavenumber; inversion to wavelength reverses it. Absorption evidence also records
whether absorption was detected, the sampled range and the source's detection
qualification. A finite-window non-detection is not an unconditional bound on the
lowest level. A missing limit is `None`, with its reason in an immutable
`missing_reasons` mapping keyed by the limit name.

Uncertainty is optional structured metadata with its value in SI, level and
provenance, separate from the kind and the modelling qualifications. Source rounding,
unknown uncertainty, a band-origin approximation and an experimental bound must remain
distinguishable. `band_edge` is a nominal derived number without an accuracy margin;
its input uncertainties and qualifications remain accessible on the limits. It uses
`e2` for `"linear"` (one common linear polarisation for all legs) and `e2_any` for
`"any"` (no polarisation restriction), with no fallback between them. Other strings
raise an `OptmatError` subclass. An eligible molecular band origin still carries the
ground-vibrational-state approximation; the result does not bound hot-band lines.

For a fitted extent, `witnesses` is an immutable tuple of point records. Each identifies
the dataset and row, process, fundamental vacuum wavelength in m, any second
wavelength (probe or Stokes), and all four signed legs in m⁻¹. Keep all ties in stable
dataset/row order; a pump wavelength alone cannot reproduce a pump-probe or CARS
quartet. The same point representation serves `max_w2_witnesses`. Only points
actually admitted to the fit count, after its ranges and exclusions; validators,
windowed-out values and scale priors do not enlarge these extents. The three maxima
summarise sampled data, not an interpolation domain or an accuracy guarantee.

A sketch (fenced as `text` until it runs; the values are argon's from the analysis):

```text
>>> from optmat.gas.gases import Gas
>>> from optmat.gas.nonlinear import hyperpolarizability_models
>>> support = hyperpolarizability_models()[Gas.Ar].support
>>> round(1e9 / support.leg_max.wavenumber, 1), support.leg_max.witnesses[0].process
(229.0, 'eshg')
>>> support.leg_max.datasets
('shelton_1990/table3_Ar_literature',)
>>> round(1e9 / support.leg_max_in_w2, 1)      # derived, sqrt(3 W2,max / 4), nm
215.9
>>> round(1e9 / support.e1.wavenumber, 1), support.e1.datasets
(106.7, ('nist_asd_2026/level_Ar_odd_J1',))
>>> round(1e9 / support.band_edge(), 1), round(1e9 / support.band_edge("any"), 1)
(189.4, 192.1)
>>> # 192.1 nm is the general-polarisation bound, not an SPM resonance.
>>> hyperpolarizability_models()[Gas.SF6].support.band_edge() is None  # E2 not curated
True
```

**Composition.** The record describes; it changes no evaluation and carries no rule.
`gamma_au` keeps its warning beyond `W2_MARGIN` times the fitted $W_2$. A propagation
code decides what to do with the facts and records the record's digest with its
output.

### 7. Implementation strategies

**Where the numbers live.**
- *(A) Computed at load* from the configurations and the curated datasets. Always
  current, but loading every dataset of a gas costs seconds, and the selections (which
  level is $E_1$) are code.
- *(B) A curated selection file and a stored derivation* (recommended):
  `data/gas/nonlinear_support.toml` names, per gas, the datasets of each limit with the
  reason, as `ionisation.toml` does; `data/gas/models/nonlinear_support/<gas>.json`
  holds the derived numbers, written by a tool and checked by a freshness test that
  re-derives them. Fast to load and reviewed as data. The selection specifies the
  candidate levels and the minimum or band-origin sum used, including which candidate
  wins; `e2_any` is the minimum of all eligible candidates, not automatically the
  $J=1$ level. The derivation identifies the fit configuration, selected dataset
  contents and schema/derivation version. Its digest changes when those inputs or
  their qualifications change, even if rounded displayed values do not.

  Reuse the analysis tool's process conversion, but do not copy its permissive
  `_values` exception handling: the production derivation must fail on an unreadable
  or unsupported fitted dataset. Recomputed maxima and the contributing point count
  must agree with the stored fit; any discrepancy is investigated rather than hidden
  by falling back to the fit's single stored maximum. Missing scientific evidence is
  explicit metadata; a broken derivation is an error.
- *(C) Inside the stored fit* (`models/nonlinear/<gas>.json`). Rejected: the stored fit
  enters the model's digest whole, so every model's digest and golden value would change
  for a description that changes no evaluation.

**Identity.**
`support` is not a field of the model and leaves its digest unchanged; the record is
itself `Identified` with its own digest, which solprop records alongside the model
digest in its manifest and resume key. Its `model_digest` must equal the requesting
model's digest, so a gas name alone cannot attach the default fit's support to a
different polynomial. A directly constructed or changed model without a matching
curated record returns `None`; this OEP adds no override input to
`HyperpolarizabilityModel`. A missing or mismatched support file for a shipped curated
model is a packaging/freshness error and raises an `OptmatError` subclass. Loading the
support must not recursively construct the same model or hash the support into the
model. Use a separate derivation version, without changing the model's evaluation
version for this metadata-only addition. An evaluator applying a cap inside optmat
was considered and is not proposed: how to use the facts is the propagation code's
decision (TC-868).

**Vocabulary and data.** No change: the quantities and the curated levels exist
(TC-862, TC-863). The polyatomic limits that are not curated are recorded as missing, not guessed.

**Cost.** One small JSON per gas; no new dependency.

### 8. Verification and validation

- **The derived leg bound:** no fitted value of any gas, and no quartet within the
  fitted $W_2$, has a leg beyond $\sqrt{3 W_{2,\max}/4}$ (the bound of section 4, tested
  on the stored data, on random quartets and on the third harmonic at $W_{2,\max}$).
- **The limits reproduce their sources:** each limit reproduces its datasets' printed
  values after the documented SI conversion or band-origin sum (NIST ASD levels;
  Huber & Herzberg's Te and 0-0 bands, with the chains of bands to the ground level).
  Cross-check helium against solprop's SEP 12, section 4.4 (two-photon 1s2s ¹S at
  a 120.3 nm fundamental, three-photon 1s2p ¹P at 175 nm, three-photon ionisation below
  151 nm); xenon's 6p[1/2]₀ at 249.6 nm against Lehmberg et al. [2].
- **Monotonicity:** $A > 0$, $B \ge 0$, $B + B' > 0$ for every stored fit, as a test that
  fails if a refit breaks it; also require $\gamma_0>0$. This tests a polynomial
  property, not a bound on physical gamma.
- **Evidence, not a property of the record:** the polynomial against every calculated
  value beyond the fitted data stays on the accuracy page (TC-867).
- **Extents and provenance:** reconstruct every fitted quartet from its row and all
  wavelength inputs. Test W2, W4, all six distinct pairs, the three maxima and their
  witnesses against an independent enumeration. Include static, mixed-sign,
  pump-probe and CARS examples, exclusions, ties and missing/unsupported inputs;
  conversion failures must fail the derivation. Assert the m⁻¹/cm⁻¹ and squared-unit
  conversions explicitly. Numeric comparisons use a stated floating-point tolerance
  (for example 1e-12 relative for these short sums), not a physical uncertainty.
- **Limits and derived bound:** test `e2_any <= e2` where both exist, candidate
  selection and chains to the ground state, unknown uncertainties and missing reasons.
  Keep absorption evidence out of `band_edge`; a single known limit must not fill
  a missing one. Test both valid polarisation values and invalid input. Verify
  $|s_a+s_b|\le2K$ including outputs, and the atomic equal-frequency $J=1$
  cancellation of [9]; do not identify the general-polarisation bound with an SPM
  pole. Rounded wavelength examples use half their last printed digit as tolerance;
  molecular approximations remain separate metadata.
- **Freshness and identity:** the stored derivation equals a fresh one; the record's
  digest is golden-tested. Tests must show that changing a witness, evidence kind,
  missing reason, uncertainty, selection or contributing dataset changes that digest;
  changing support alone leaves the model's digest and evaluated gamma unchanged.
  Test nested immutability, unmatched user models, corrupted/missing shipped records
  and the contents of an installed wheel/sdist. The full gate and release test are
  required for implementation; no refit is needed while all fit inputs stay unchanged.
- **Parity with solprop:** SEP 6's Phase 1 reads the record in place of its own copies of
  the numbers.

### 9. Documentation

A subsection of {doc}`../../nonlinear` on the record, with the doctest of section 6 and a
link from the analysis there; the record on {doc}`../../identity`; the public-API
reference; a line on {doc}`../../propagation`. The existing support and accuracy text was
reconciled with section 4 after the review (TC-869): a sufficient bound is not an
inevitable resonance, finite-window absorption evidence is not a level, and the
calculated comparisons are not a universal lower-bound proof. Phase 1 keeps the public
record's labels consistent with that. A gallery example is not needed.

### 10. Compatibility and impact

- No recommended value, default, stored fit or evaluated gamma changes.
- No model digest changes (`support` is not a field of the model; see *Identity*); the
  public-API snapshot gains the new names.
- New data files and a freshness test; the public-surface snapshot and the release test
  change.
- No solprop parity fixture changes until SEP 6 uses the record.

### 11. Alternatives considered

- **A record with a recommended rule** (a per-leg cap, a refusal at the band edge, a
  status per gas), as first drafted. Withdrawn: optmat states facts, and their use is a
  modelling decision that belongs to the propagation code (TC-868). The comparison of
  caps that motivated it stays on the accuracy page as evidence (TC-867).
- **An evaluator of capped quartets in optmat.** Not proposed, for the same reason; it
  would also change every model's digest.
- **Publish the numbers only in the manual.** solprop would copy them, and a refit or a
  new source would silently go stale on its side.
- **Warn only.** A warning is not a fact a propagation code can build a policy on.

### 12. Phases

1. **The record.** Its data are curated (TC-862, TC-863) and its numbers are computed
   by `tools/fit/nonlinear_support.py` (TC-864 to TC-868). Deliverables: the selection
   file and the stored derivation, `HyperpolarizabilitySupport` and
   `HyperpolarizabilityModel.support`, the freshness and golden tests, the manual
   subsection and corrections in section 9, the register entries and changelog entry.
   Verifies section 8 except solprop's parity, which belongs to SEP 6 and does not
   block optmat's phase. Its report records any unverified source or deferred limit.

### 13. Open questions

solprop's SEP 6 review answered the first draft's questions (2026-10-07):

- **Q1-Q3 (the use of the facts)** are solprop's: a lower bound below the band edge and
  best effort beyond it; the xenon and CO mixed quartets noted only; solprop checks the
  region itself, so no band-check function is needed. Closed for optmat's ownership
  decision; section 4 qualifies the evidence and does not certify solprop's
  lower-bound interpretation.
- **Q4, Q5 (the polyatomic molecules, O3):** the record states the facts as they are: no
  two-photon state is curated for CO2, CH4, SF6 or N2O, and O3 absorbs throughout its
  range. Closed.
- **Q6 (O2's spin-forbidden a ¹Δg and b ¹Σg⁺ states):** a material judgement, optmat's:
  not counted among the limits (TC-865, TC-868); solprop agrees. Closed.
- **Q7 (units):** m⁻¹ (SI). Closed.
- **Q8 (an accuracy margin below each resonance):** published as information if optmat
  has evidence for one, never as a rule. Deferred to the roadmap by Phase 1: optmat
  has none today; this feature supplies no accuracy margin.
- **Q9 (Orr & Ward's original):** resolved by removing the citation, whose full text was
  not checked. Shelton & Rice's printed eq. 14 (p. 8) carries the whole denominator
  structure section 4 uses, with the permutations over all four frequencies, and is
  checked (the author, 2026-10-07).
- **Q10 (existing analysis qualifications):** the manual, the analysis tool's labels
  and the register were reconciled with sections 4-8 after the review (TC-869). Phase 1
  exports no unqualified absorption endpoint, fallback band edge or universal
  lower-bound wording as a public fact.

### 14. References

1. D. P. Shelton and J. E. Rice, "Measurements and calculations of the
   hyperpolarizabilities of atoms and small molecules in the gas phase," *Chem. Rev.*
   **94**, 3 (1994), eq. 14 on p. 8 and eq. 17 on p. 10. doi:10.1021/cr00025a001 (curated,
   `shelton_1994`).
2. R. H. Lehmberg, C. J. Pawley, A. V. Deniz, M. Klapisch and Y. Leng, "Two-photon
   resonantly-enhanced negative nonlinear refractive index in xenon at 248 nm," *Opt.
   Commun.* **121**, 78 (1995). doi:10.1016/0030-4018(95)00433-9 (curated,
   `lehmberg_1995`).
3. K. P. Huber and G. Herzberg, *Molecular Spectra and Molecular Structure IV. Constants
   of Diatomic Molecules* (Springer, 1979). doi:10.1007/978-1-4757-0961-2 (not curated;
   read through the NIST Chemistry WebBook; checked on Crossref).
4. A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team, NIST Atomic Spectra Database
   (version 5.12). doi:10.18434/T4W30F (curated, `nist_asd_2026`).
5. D. M. Bishop and J. Pipin, *J. Chem. Phys.* **91**, 3549 (1989).
   doi:10.1063/1.456886 (curated, `bishop_1989`).
6. D. Ityaksov, H. Linnartz and W. Ubachs, *Mol. Phys.* **106**, 2471 (2008).
   doi:10.1080/00268970802570334 (curated, `ityaksov_2008`).
7. W. F. Chan, G. Cooper and C. E. Brion, *Chem. Phys.* **180**, 77 (1994).
   doi:10.1016/0301-0104(93)E0386-A (curated, `chan_1994`).
8. A. Serdyuchenko, V. Gorshelev, M. Weber, W. Chehade and J. P. Burrows, *Atmos. Meas.
   Tech.* **7**, 625 (2014). doi:10.5194/amt-7-625-2014 (curated, `serdyuchenko_2014`).
9. D. DeMille, D. Budker, N. Derr and E. Deveney, "Search for Exchange-Antisymmetric
    Two-Photon States," *Phys. Rev. Lett.* **83**, 3978 (1999).
    doi:10.1103/PhysRevLett.83.3978 (not curated; authors' manuscript
    [arXiv:physics/9906025v2](https://arxiv.org/abs/physics/9906025v2), eqs. 1-4 checked
    in the review; the author accepts that check, 2026-10-07).

## Part II — Implementation record

### Phase 1 — The record (2026-10-07)

**Delivered.**
- `optmat.gas.nonlinear`: `HyperpolarizabilitySupport` (`Identified`, kind
  `gas_hyperpolarizability_support`), `SupportLimit`, `SupportPoint`,
  `SupportUncertainty`, `SupportCandidate`, `AbsorptionEvidence`,
  `SUPPORT_LIMIT_KINDS`, `SUPPORT_DERIVATION_VERSION` (1), and the read-only property
  `HyperpolarizabilityModel.support`. The stored records are read once and cached.
- `optmat.gas._nonlinear_support` (private): the derivation, `signed_legs`,
  `fitted_points` and `derive_support`; `tools/fit/nonlinear_support_record.py` writes
  the records.
- `data/gas/nonlinear_support.toml`: the selection, with each limit's candidates,
  register entries and qualifications, and the reason for each missing limit.
- `data/gas/models/nonlinear_support/<gas>.json`: the 15 records.
- Tests: `tests/test_nonlinear_support_record.py` (237 tests); the golden file gains a
  record per gas (its numbers and digest); the public-API snapshot gains the new names;
  the snippets of *Third-order nonlinearity* and *Using optmat from a propagation code*
  read the record, the latter also from the installed wheel and sdist.
- The manual: *The support record* on {doc}`../../nonlinear` (and the command after a
  refit); *Records about a model* on {doc}`../../identity`; the third order under *Where a
  model holds* and the surface table on {doc}`../../propagation`. `AGENTS.md`: the
  re-derivation after a hyperpolarizability refit.

**Deviations, and why.**
1. Section 6 names `HyperpolarizabilitySupport` and `SupportLimit`; the record's nested
   values are public frozen dataclasses too (`SupportPoint`, `SupportUncertainty`,
   `SupportCandidate`, `AbsorptionEvidence`), since a caller receives them. The
   derivation and the JSON reader and writer stay private.
2. Beyond section 6's table, the record holds `fitted_points` (the count the derivation
   checks against the stored fit), `derivation_version`, `input_digests` (the
   selection's, the fit configuration's and each selected dataset's SHA-256, which
   section 7 requires to enter the digest) and, on a limit taken as the smallest of
   several, its `candidates` (section 7: "including which candidate wins").
3. An unknown polarisation raises `DataError`, the `OptmatError` subclass that
   `frequency_sums` raises for an unknown process.
4. The ionisation energy is a limit of kind `level`, the kind of a threshold energy
   taken as its source gives it; its own kind (adiabatic, vertical) is a qualification.
5. A fitted extent and absorption evidence carry no uncertainty (None, not
   applicable); levels, band origins and ionisation energies always carry one, stated
   or explicitly unstated.
6. Absorption evidence's `sampled` window is the curated dataset's span, which can be
   narrower than the measurement (O3's table ends at 1046.8 nm, the paper's range at
   1100 nm); the qualification says so.
7. The record's identity status is `validated`: it is curated and tested. It is the
   record's standing as data, not a judgement on any quartet.
8. The analysis tool `tools/fit/nonlinear_support.py` keeps its own selection for the
   manual's tables; a test requires every number of the record to agree with it (Q11).
9. The manual's remaining citation of Orr & Ward, removed from Part I in review (Q9),
   is replaced by Shelton & Rice's eq. 14.

**Decisions.** TC-870: the selection (TC-865 as qualified by TC-869), the extents and
their witnesses (ties within 1e-12 relative), and the agreement with the stored fit
required exactly. No recommended value, stored fit or evaluated gamma changed.

**Register.** TC-870.

**Verification.**
- `uv run poe check`: all eight steps passed (18,201 tests, one expected skip, 100 %
  statement and branch coverage, the documentation built).
- `uv run poe test-release`: the wheel ships the selection and the records, and the
  propagation page's snippet reads `support` from the installed wheel and sdist.
- The records reproduce the analysis of the manual: extents and point counts to
  1e-12, every limit to 1e-9 (the tool's 10-digit hc/e), every band edge to its printed
  0.1 nm. Levels are their printed NIST values times 100; band origins the sums of
  their printed 0-0 bands; helium agrees with SEP 12, section 4.4, and xenon's
  6p[1/2]0 with Lehmberg et al.; every ionisation energy is `optmat.gas.ionisation`'s.
- Every model digest and golden value is unchanged; any change to a witness, kind,
  missing reason, uncertainty, selection digest or dataset changes the record's digest.
- A user's model has no record; a curated model without its record, a corrupted or
  misnamed record, an older derivation version, an unreadable fitted dataset, and a
  count or largest W2 that disagrees with the stored fit each raise `DataError`.

**Sources checked and not checked.** No new source. The limits come from TC-862,
TC-863 and the curated absorption datasets, checked there. The molecular one- and
two-photon assignments remain the curators' reading of the symmetries (TC-865).

**Limitations.** CO2 has no E1, and no polyatomic molecule has an E2 (TC-865, open):
their band edges are None. N2O's support is one pump-probe value. The record gives no
linewidths, tensor amplitudes or accuracy margin.

**Open questions.**
- Q11: should the analysis tool read the record's selection, leaving one copy?
  Today a test keeps the two equal.
- Q12: Final status, and the move to the archive, are the maintainers' approval.
- Q8 is on the roadmap (an accuracy margin below each resonance).

**Handoff.** solprop's SEP 6 Phase 1 reads `HyperpolarizabilityModel.support` in place
of its own copies, and stores `support.digest` beside the model's. After a refit or a
change to the selection or a selected dataset, re-run
`tools/fit/nonlinear_support_record.py` for the gas.

## Part III — History

- 2026-10-07: Drafted, for solprop's SEP 6 Q2, alongside the data change that curates
  the resonance levels and adds the analysis to the manual (TC-862 to TC-867).
- 2026-10-07: Revised to cite that change, now on `main`; submitted for review.
- 2026-10-07: Revised to facts only, after solprop's SEP 6 review and the maintainer's
  decision (TC-868): the recommended rule, the band check and the status are withdrawn;
  the two derived quantities stay, labelled as derived; the open questions record
  solprop's answers.
- 2026-10-07: Part I reviewed and revised by GPT-6 (AI). John Travers confirmed his
  completed human review; review status set to Reviewed, lifecycle status still Draft.
  Corrected the general-polarisation band interpretation and qualified absorption
  endpoints, uncertainty and the claimed physical lower bound. Specified row-level
  witnesses, units, missing evidence, model binding, digest dependencies and failure
  behaviour. Reviewed the current model/identity code, analysis tool, source records,
  stored NIST evidence and manual tables; checked Shelton & Rice's printed eqs. 14 and
  17, Ityaksov et al.'s stated non-detection, Serdyuchenko et al.'s measurement coverage
  and DeMille et al.'s two-photon cancellation. No new data were transcribed and no fit
  was changed. The original Orr & Ward paper remains unchecked (Q9); the molecular
  selection assignments remain the curators' assumptions, not independently verified
  transition-strength calculations. Q8 and the corrections to the existing analysis
  have explicit destinations in Phase 1. Reviewed does not imply accepted.
- 2026-10-07: Review validation: `git diff --check`, codespell on this proposal,
  `tests/test_proposals.py`, `tests/test_register.py`,
  `tests/test_nonlinear_support_page.py`, `tests/test_nist_asd.py` and
  `tests/test_nist_webbook_states.py` passed (133 tests, one expected Process-OEP
  skip). Confirmed positive gamma0 for all 15 stored models. `uv run poe docs` passed
  with all 53 gallery examples executed. This is a proposal-only change: no full
  package gate, extraction rerun, refit or release test was required or run.
- 2026-10-07: The author accepted the review revisions without a further review.
  Orr & Ward's unchecked citation removed in favour of Shelton & Rice's eq. 14 (Q9);
  DeMille et al.'s manuscript check accepted; the manual, tool and register reconciled
  (TC-869, Q10); section 10's reason for the unchanged model digests corrected.
- 2026-10-07: Accepted by the author, John Travers. Part I is frozen.
- 2026-10-07: In progress. Phase 1 implemented on `feature/oep-0005-support-record`
  and reported in Part II (TC-870): `uv run poe check` passed (18,201 tests, one
  expected skip, 100 % statement and branch coverage, the documentation built), and
  `uv run poe test-release` passed (4 tests). No model digest or golden value changed.
- 2026-10-07: Phase 1 merged to `main` (337f811). Final, approved by the author, John
  Travers. Q8 and Q11 are on the roadmap; Q12 is this approval.
