# OEP 4 — The refractive index of gas mixtures at finite density

| | |
|---|---|
| **OEP** | 4 |
| **Title** | The refractive index of gas mixtures at finite density |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 |
| **Status** | Final |
| **Review status** | Reviewed |
| **Reviewers** | GPT-6 (AI); John Travers (human) |
| **Type** | Feature |
| **Created** | 2026-09-30 |
| **Resolution** | 2026-09-30: accepted by the author, John Travers, after the human and AI reviews; merged directly to `main` |
| **Requires** | [OEP 1](oep-0001-propagation-interface.md): the equation-of-state object, identities and the public surface |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | Opus 5.5 (AI), for John Travers. Phases 1, 2 and 3 merged from `feature/oep-0004-optics`, `feature/oep-0004-density` and `feature/oep-0004-example` |
| **Manual pages** | {doc}`../../mixtures`; {doc}`../../conventions`; {doc}`../../data`; {doc}`../../auto_gallery/accuracy/plot_mixture_he_ar` |

## Part I — Proposal

### 1. Abstract

optmat gives the index of a pure gas at any density: Lorentz–Lorenz with the second
refractivity virial coefficient, at the real-gas density of its equation of state. It has
no law for a mixture. This OEP proposes one, in two parts that belong together:

- **the optics of a mixture:** the Lorentz–Lorenz function summed over the species and
  their pairs, with the unlike-pair refractivity virial coefficients, inverted once;
- **the density of a mixture:** a mixture equation of state, which carries the
  unlike-pair interactions that the sum of the constituents' densities at their partial
  pressures leaves out.

The second matters more. For the mixtures a gas-filled fibre uses (He/Ar, Ne/Ar, air), the
partial-pressure rule's density differs from the multi-fluid mixture models' by 3 × 10⁻⁴
to 1.6 × 10⁻³ at 10 bar, and the rule does not preserve the declared composition. Any
plausible unlike-pair refractivity term changes $n - 1$ by less than 10⁻⁴ (§3). These are
differences between models. Near ambient the multi-fluid models deviate from the
measured He–Ar and Ne–Ar densities by several tenths of a per cent, more than those
differences. So Phase 2 compares every density method with measured mixture densities
and cross second virial coefficients, and chooses each mixture's default from that. The work is what
solprop's SEP 11, the finite-density gas law, needs from optmat.

Success is judged in three ways:

- the law reduces exactly to optmat's pure-gas law for one species;
- the mixture density agrees with CIPM-2007's air formula near ambient conditions within
  that formula's stated uncertainty, and with measured cross second virial coefficients;
- every unlike-pair term carries a source or a declared estimate with its uncertainty,
  and its effect on $n - 1$ is reported over a stated domain.

### 2. Motivation and scope

A gas-filled fibre's zero-dispersion and phase-matching pressures are often set with a
mixture, helium with argon or neon with argon, at several to tens of bar. solprop
propagates labelled populations of species along a fibre, each with its own density. It
now uses the dilute law, $\chi = \sum_s N_s\alpha_s/\varepsilon_0$, and the
partial-pressure density rule. SEP 11 proposes to replace the dilute law with optmat's
finite-density law, and asks optmat for the mixture law with its unlike-pair terms and
their domains (SEP 11, §6, *Upstream*). OEP 1 left this out of its scope (its §2).

**In scope:**

- the Lorentz–Lorenz function of a mixture in the species' densities, with the like- and
  unlike-pair refractivity virial terms, its inversion, and its derivatives in frequency
  and in each density;
- a record for unlike-pair coefficients, curated or declared, with their status and
  uncertainty;
- the density of a mixture from a mixture equation of state, and the inverse;
- the domain of each, and the sensitivity of $n - 1$ to the terms that are not known.

**Out of scope:**

- the third refractivity virial coefficient, neglected as for the pure gases now (§4.1);
- density effects on the nonlinear response (local-field factors, pair terms of
  $\chi^{(3)}$), on Raman lines or on ionisation;
- solprop's consumers of the law (its caches, fast paths and expansions): SEP 11's;
- the like-pair coefficients that optmat has not curated (O2, D2, N2O, CO, O3), which
  are roadmap items; this OEP records their omission and its effect.

### 3. Context

optmat's pure-gas law is in `optmat.gas.lorentz_lorenz` and `optmat.gas.virial`
({doc}`../../conventions`, *Polarizability and refractive index*):

$$
\frac{n^2-1}{n^2+2} = A_R\rho + B_R\rho^2, \qquad A_R = \frac{N_A\alpha}{3\varepsilon_0},
$$

with $\rho$ the molar density and $B_R$ chosen per gas in
`data/gas/refractivity_virial.toml` (TC-170, TC-178): Achtermann et al. at
633 nm for He, Ar, H2, N2, CH4, CO2 and SF6; Orcutt and Cole's static dielectric
coefficient for Ne and Kr; Burns et al. for Xe. O2, D2, N2O, CO and O3 have none. The
code writes it as $y = N\alpha(1 + b\rho)/(3\varepsilon_0)$ with $b = B_R/A_R$ constant,
so that the pair term disperses like the linear one.

The density is `optmat.gas.eos.number_density`, a pure gas's real-gas density from its
CoolProp equation of state; OEP 1 (§4.2, §6.5) makes it an object with an identity, a
phase domain and an inverse, and adds CoolProp's pseudo-pure dry air for `Gas.Air`.

**Unlike-pair data curated today.** Three sources give static dielectric data of
mixtures, recorded in notes and not extracted:

- Bose and Cole, CO2–Ar (`bose_1970`): the cross dielectric coefficient
  $B_{\varepsilon,12} = 14.3$ cm⁶/mol² at 322.8 K, of which the induced part, the one
  that survives at optical frequencies, is "assumed" (2.4) rather than measured; and the
  measured cross pressure virial coefficient $B_{p,12} = -30.1$ cm³/mol;
- Hourri et al., CO2–SF6 (`hourri_1997`): the optical cross term $B_{R,12} = 6.57$
  cm⁶/mol², which is the geometric mean of the pure $B_R$, not a measurement;
- Copeland and Cole, CF3H–Ar and CFH3–Ar (`copeland_1976`), species optmat does not
  carry.

So no measured optical unlike-pair coefficient of any pair optmat carries is in hand.

**The sizes involved**, computed for this OEP at `ea0dcaf` at 800 nm and 294 K, with
CoolProp 8.0.0's multi-fluid mixture models as the reference density:

| Mixture | Pressure | Partial-pressure density, relative difference | Lorentz–Lorenz against the dilute law, in $n - 1$ | Unlike refractivity term, 0 against the like-pair mean, in $n - 1$ |
|---|---:|---:|---:|---:|
| He/Ar 50/50 | 10 bar | +1.6 × 10⁻³ | +1.0 × 10⁻³ | +7 × 10⁻⁵ |
| He/Ar 50/50 | 30 bar | +5.5 × 10⁻³ | +3.1 × 10⁻³ | +2.2 × 10⁻⁴ |
| Ne/Ar 50/50 | 10 bar | +3.4 × 10⁻⁴ | +1.1 × 10⁻³ | +6 × 10⁻⁵ |
| Ne/Ar 50/50 | 30 bar | +1.7 × 10⁻³ | +3.4 × 10⁻³ | +1.7 × 10⁻⁴ |
| air (N2/O2/Ar) | 10 bar | −1.6 × 10⁻³ | +1.8 × 10⁻³ | +2 × 10⁻⁵ |
| air (N2/O2/Ar) | 30 bar | −4.1 × 10⁻³ | +5.5 × 10⁻³ | +5 × 10⁻⁵ |

The columns are relative changes, of the total density or of $n - 1$, not absolute
increments. The reference for the density column is CoolProp's multi-fluid model at the
same pressure, temperature and fixed composition. The unlike column compares a zero
unlike term with the mean of the like terms' $B$ (§7.2), with the like terms kept, at the
multi-fluid density. Air is the CIPM major fractions (0.780848, 0.209390, 0.009332),
renormalised to one. These are differences between models, not errors established
against experiment. Near ambient, Tkaczuk et al.'s models [4] deviate from the measured
mixture densities on average by (TC-742):
- 0.08–0.26 % for He–Ne (Kestin et al. 1981 at 298 K; Vogl and Hall 1972);
- 0.72 % for He–Ar (Blancett et al., 223–323 K);
- 0.43 % for Ne–Ar (Kestin et al. 1981 at 298 K).

These are more than the models' differences from the partial-pressure rule in the
table. Which method is closer to experiment is Phase 2's comparison (§7.4).

The density difference is as large as the whole Lorentz–Lorenz correction, and ten to
thirty times the uncertainty of the unknown unlike-pair optics.

The partial-pressure rule also changes the composition. With each species' own equation
of state at $x_s p$, the number fractions differ from $x_s$. For 50/50 He/Ar at 294 K the
helium fraction becomes 0.49856 at 10 bar and 0.49574 at 30 bar. Against the multi-fluid
density at fixed composition, $n - 1$ at 800 nm then moves by 3.8 × 10⁻³ and
1.2 × 10⁻² (GPT-6's review, R3). The fallback of §4.2 preserves the composition. For air, CoolProp's pseudo-pure
fluid and its multi-fluid N2/O2/Ar mixture differ by 2.7 × 10⁻⁵ at 1 bar and 2.0 × 10⁻⁴
at 10 bar.

The roadmap's *solprop* section lists this law as an informal item; §12 of this OEP
covers it.

### 4. Theory

#### 4.1 The optics of a mixture

Write $n_s$ for the number density of species $s$ and $\rho_s = n_s/N_A$ for its molar
density. In the low-density cluster expansion, the Lorentz–Lorenz function of a gas is the
sum of the isolated molecules' polarizabilities and the excess polarizability of each
interacting pair, averaged over the pair's separation [1]:

$$
y \equiv \frac{n^2-1}{n^2+2}
= \frac{1}{3\varepsilon_0}\left[\sum_s n_s\alpha_s
+ \frac12\sum_{s,t} n_s n_t \int \Delta\alpha_{st}(R)\,e^{-U_{st}(R)/k_BT}\,d^3R\right]
+ O(n^3),
$$

where $\Delta\alpha_{st}(R)$ is the isotropic part of the polarizability of an $s$–$t$ pair
at separation $R$ minus $\alpha_s + \alpha_t$, and $U_{st}$ is the pair potential. The
factor ½ counts each unordered pair once in the double sum over ordered species; for
$s \ne t$ both orders appear, giving $n_s n_t$ pairs, and for $s = t$ there are
$n_s^2/2$. The integral is classical, as Buckingham and Pople's is (TC-738). Quantum
corrections are below 1 % near room temperature for helium and neon (0.9 % and 0.5 % at
300 K [5]), grow at low temperature, and are included in the measured coefficients. For
molecules the integral also averages over relative orientations. In molar terms,

$$
y = \sum_s A_s\rho_s + \sum_{s,t} B_{st}\,\rho_s\rho_t,
\qquad
A_s = \frac{N_A\alpha_s}{3\varepsilon_0},
\qquad
B_{st} = B_{ts} = \frac{N_A^2}{6\varepsilon_0}\int \Delta\alpha_{st}\,e^{-U_{st}/k_BT}\,d^3R.
$$

For one species $B_{ss} = B_R$, the pure gas's coefficient, so the mixture law reduces to
optmat's pure law exactly. For mole fractions $x_s$ at total molar density $\rho$,
$B_\text{mix} = \sum_{s,t} x_s x_t B_{st}$, the familiar quadratic mixing rule, which
follows from the pair count and is not an assumption.

**Three consequences.**

- **Invert once.** $\chi = 3y/(1-y)$ and $n - 1 = \chi/(\sqrt{1+\chi} + 1)$ are applied to
  the whole mixture's $y$. Summing the species' indices, or their $\chi$ after inverting
  each, is wrong at the same order as the pair terms.
- **Labelled populations combine first.** Two labelled populations of one species are one
  species in the pair sum; otherwise splitting a population into two would change its
  refractivity.
- **The unlike-pair terms are data.** $B_{st}$ depends on the $s$–$t$ pair polarizability
  and potential, and cannot be deduced from $B_{ss}$ and $B_{tt}$. A combining rule is an
  estimate. Hourri et al. used the geometric mean for CO2–SF6 (§3), which is undefined
  when the like terms differ in sign, as helium's ($-0.068$ cm⁶/mol²) and argon's (1.76)
  do.

**Dispersion of the pair terms.** optmat's pure law makes the pair term disperse like the
linear one, $B_{ss}(\omega) = b_{ss}A_s(\omega)$. The mixture law keeps that convention,
with $B_{st}(\omega) = b_{st}\,\tfrac12[A_s(\omega) + A_t(\omega)]$ for an unlike pair. It
is an approximation to a small term: Burns et al.'s $B_R$ of xenon at four lines spreads
from 14.9 to 25.5 cm⁶/mol² (TC-178), so the dispersion of the pair term is not resolved by
the data either way. For He, Ne and Ar, Garberoglio and Harvey's calculation [5] puts
the pair term's dispersion, from static to 633 nm, at a few per cent. The estimate of
§7.2, the mean of the like terms, disperses as
$\tfrac12[B_{ss}(\omega) + B_{tt}(\omega)]$, each like term with its own $A$.

**The neglected third term.** The expansion is truncated after the pair term, as for the
pure gases now. The domain of §4.3 bounds its omission with measured third coefficients.

**Derivatives.** With $\partial\chi/\partial y = 3/(1-y)^2$:

$$
\partial_\omega\chi = \frac{3}{(1-y)^2}\left[\sum_s\rho_s\,\partial_\omega A_s
+ \sum_{s,t}\rho_s\rho_t\,\partial_\omega B_{st}\right],
\qquad
\frac{\partial\chi}{\partial\rho_s} = \frac{3}{(1-y)^2}\left[A_s + 2\sum_t B_{st}\rho_t\right].
$$

$\partial_\omega A_s$ is OEP 1's `dpolarizability_domega` times $N_A/3\varepsilon_0$. Both
derivatives are closed forms.

**Smoothness in density.** $y$ is a quadratic in the densities and $\chi$ has its only
singularity at $y = 1$, twenty to two hundred times the largest density of interest. A
propagation code that expands its linear operator in density loses nothing: the
Chebyshev degree needed to represent $\sqrt{1+\chi} - 1$ to 10⁻¹³ over zero to full
density is 5 for argon at 30 bar, 6 for N2 at 100 bar and 7 for xenon at 60 bar with
Lorentz–Lorenz, against 5, 7 and 8 with the dilute law (computed for this OEP at 800 nm).

#### 4.2 The density of a mixture

To second order in density, a mixture's compressibility factor is

$$
Z = 1 + \rho\sum_{s,t}x_s x_t B_{st}^{(p)}(T) + O(\rho^2),
$$

with $B^{(p)}_{st}$ the second (pressure) virial coefficients, like and unlike, and the same
pair count as §4.1 (Gernert's Eq. 5.11 [9]). A multi-fluid model's second virial is not
quadratic in composition in general: its effective cross coefficient depends on
composition [9] (TC-743). The partial-pressure rule puts each species at $p_s = x_s p$ with
its own equation of state, which to the same order gives
$Z_\text{pp} = 1 + \rho\sum_s x_s^2 B^{(p)}_{ss}$. It is therefore the mixture law with
every unlike-pair coefficient set to zero, and at fixed pressure and temperature

$$
\frac{\rho}{\rho_\text{pp}} = 1 - \rho\sum_{s \ne t}x_s x_t B^{(p)}_{st} + O(\rho^2).
$$

An unlike-pair pressure coefficient of 10 cm³/mol (CO2–Ar's is −30.1 cm³/mol,
`bose_1970`) makes the term 2 × 10⁻³ for an equimolar mixture at 10 bar
($\rho \approx 400$ mol/m³): the size in §3's table. Beyond first order the rule also
changes the composition (§3).

**The fallback.** Where no qualified multi-fluid model covers every pair, the pressure is

$$
p(\rho, T) = \sum_s p_s(x_s\rho, T) + RT\rho^2\sum_{s \ne t}x_s x_t B^{(p)}_{st}(T),
$$

with $p_s$ species $s$'s own equation of state at its partial molar density $x_s\rho$, and
$B^{(p)}_{st}$ the curated cross second virial coefficients.
- It preserves the composition: $N_s = x_s N$ exactly.
- It reduces to the pure equation of state at $x_s = 1$, and to
  $Z = 1 + \rho\sum_{s,t}x_s x_t B^{(p)}_{st}$ at low density.
- The pressure is explicit in the densities, which suits a propagation code that holds
  labelled densities. The density at a given pressure is a one-dimensional root,
  monotone in the gas phase.
- With every $B^{(p)}_{st} = 0$ it is the composition-preserving form of the
  partial-pressure rule: the method `partial-density`.
- It leaves out the unlike third virial terms. It is qualified while
  $|\rho\sum_{s\ne t}x_s x_t B^{(p)}_{st}| \le 10^{-2}$, with a truncation uncertainty
  of the square of that sum. That is a stated assumption, checked against the
  multi-fluid density on the pairs that have both (§8). Beyond it, the fallback is
  provisional.

A **multi-fluid mixture equation of state** carries these interactions, and the pair's
departure beyond second order, through binary parameters fitted to mixture data; it
reduces to each pure fluid's equation of state at $x_s = 1$. CoolProp implements the
GERG-2008 form [3] and later binary models. Among the 91 pairs of the fourteen gases
optmat carries (He, Ne, Ar, Kr, Xe, H2, D2, N2, O2, CH4, SF6, N2O, CO2, CO), CoolProp 8.0.0
has binary parameters for 39, among them He–Ne, He–Ar and Ne–Ar (Tkaczuk et al. [4]), and
N2–O2, Ar–N2 and Ar–O2. It has none for He–Kr, He–Xe, Ne–Kr, Ne–Xe, Ar–Xe or most pairs
with SF6 and N2O; there it refuses the mixture (Ar–Xe fails at the first call), and the density falls
back (§7.4).
- The 39 records cover different models and levels of evidence: 17 `Kunz-JCED-2012`,
  ten `Gernert-Thesis-2013` [9], nine `Bell-JCED-2016` and three `Tkaczuk-JPCRD-2020`
  (GPT-6's review, R6).
- Six are placeholders, with every reducing parameter one and no departure function:
  He–H2, He–O2, He–CO, Ar–H2, H2–O2 and O2–CO. These are GERG-2008's Lorentz–Berthelot
  pairs, whose parameters were not fitted to data of the pair [3] (TC-740).
- Air's pairs are GERG-2008's N2–O2, with fitted reducing functions, and Gernert's N2–Ar
  and O2–Ar [9]. Gernert states gas-phase density uncertainties below 0.1 %, 0.5 % and
  1.0 % for them. They are not Lemmon et al.'s air mixture model [6], whose pseudo-pure
  air equation is OEP 1's `Gas.Air` (TC-741).

#### 4.3 The domain

The law is second order in density, and its unlike-pair terms are known only as well as
their sources. Two domains are kept apart.

**The optical domain** is a density bound per mixture, the lowest of:

- **the third term.** The density at which the neglected third term reaches 10⁻³ of $y$.
  - Its ratio to the linear term is $C_R\rho^2/A_R$, with $C_R$ the third refractivity
    virial coefficient.
  - It uses the measured $C_R$ where they exist: Achtermann et al.'s for Ar, N2, CH4,
    CO2, SF6 and H2, and their bound for He.
  - For a mixture it takes the largest constituent $|C_R/A_R|$, since the cross third
    coefficients are unknown.
  - Where $C_R$ is not measured (Ne, Kr, Xe, O2, D2, N2O, CO, O3), a declared
    provisional scale is used: $|C_R| \le 600\,B_R^2/A_R$. It comes from the largest
    measured ratio, N2's 560. The ratio is not of order one (R1).
  - A gas with neither $B_R$ nor $C_R$ has no third-term bound, and its mixtures are
    provisional.
  - For N2 and Ar the bound is about 6700 and 7000 mol/m³, near 160 bar at 294 K.
- **a declared or unknown unlike term.** The density at which its uncertainty reaches
  10⁻³ of $n - 1$ (§7.2).

In this law the pair and third terms disperse like the linear term, so their ratios to it
do not depend on wavelength, except the estimate of §7.2. The bound is therefore taken as
the smallest over the constituents' qualified intervals, and holds at every wavelength
the mixture accepts. The coefficients are not modelled in temperature (§7.5), so the
optical bound does not depend on it. Beyond the bound the optics warn
(`ModelRangeWarning`), as the pure models do.

**Phase admissibility** belongs to the equation of state, which binds a temperature at
each state (OEP 1, §4.2).
- The multi-fluid model checks the mixture's own phase stability, with CoolProp's
  mixture flash.
- The fallback checks each pure fluid's phase domain at its partial density. That is a
  necessary condition, not a sufficient one, and the page states that its mixtures are
  not checked for phase separation.
- A state outside the phase domain raises `StateError`, as in OEP 1.

In wavelength, a mixture is valid in the intersection of its constituents' qualified
intervals (OEP 1, §4.1 and §4.4). Today each constituent has one interval.
- After [OEP 3](oep-0003-infrared-windows.md), a mixture with CO2, CH4, N2O or CO
  inherits their absorption bands.
- The intersection is then conservative for a trace constituent, whose band absorbs in
  proportion to its mole fraction. OEP 3 states when a mixture may widen it (its §13
  Q13), for air and for these mixtures alike. All three of these must pass at the
  mixture's mole fractions:
  - the absorption weighted by mole fraction;
  - no constituent pole in the wider interval;
  - each constituent's representation bound times its mole fraction within the
    tolerance.
- This law then takes that rule. The rule changes where the law warns, not its values.
  To cross a trace constituent's kept band, such as CO2's laser bands in air, needs a
  different representation of that constituent. That is a different constituent model,
  with its own identity, and changes the values.

### 5. Primary sources and data

| Needed for | Source | State |
|---|---|---|
| Like-pair refractivity virials | Achtermann et al.; Orcutt and Cole; Burns et al. | Curated (`achtermann_1991`, `orcutt_1967`, `burns_1986`) |
| Unlike-pair refractivity virials | Measured optical or low-frequency coefficients of the pairs of interest, or ab initio pair polarizabilities | None in hand (§3). Garberoglio et al.'s review of ab initio fluid properties (`garberoglio_2023`, curated) gives like-pair dielectric and refractivity virials for He, Ne and Ar and none for mixtures or molecules; Garberoglio and Harvey's path-integral calculation [5] (`garberoglio_2020`, curated, TC-739) is of He, Ne and Ar alone. To be sought in Phase 1, for He–Ar, Ne–Ar, He–Ne and N2–O2 first |
| Theory of the dielectric virial | Buckingham and Pople [1] | Curated (`buckingham_1955`, TC-738): §4.1's form checked against the paper |
| Third refractivity virials (the optical domain) | Achtermann et al., Table II | Curated (`achtermann_1991`), in notes: to be extracted in Phase 1 |
| Mixture equation of state | GERG-2008 [3], Tkaczuk et al. [4] and Gernert [9], as CoolProp implements them: He–Ar, Ne–Ar and He–Ne `Tkaczuk-JPCRD-2020`; air's N2–O2 GERG-2008's and N2–Ar, O2–Ar Gernert's. Lemmon et al. [6] for air's pseudo-pure equation | Curated (`kunz_2012`, `tkaczuk_2020`, `gernert_2013`, `lemmon_2000`; TC-740 to TC-743). Each pair's record names its model by CoolProp's reference key |
| Cross second virial coefficients (validation and fallback) | Measured values per pair, from the primary papers the models cite. For He–Ne, He–Ar and Ne–Ar: Brewer and Vaughn (1969), Martin et al. (1982), Vogl and Hall (1972), Blancett et al. (1970), Kalfoglou and Miller (1967) [4]. For N2–O2, N2–Ar and O2–Ar: those Gernert lists, among them Martin et al. and Vatter et al. (1996) [9] | Not in hand; to be requested and curated in Phase 2. The Landolt–Börnstein compilation [7] could not be obtained |
| Mixture densities near ambient (validation) | Kestin et al. (1981), Vogl and Hall (1972) and Blancett et al. (1970) for the noble-gas pairs [4]; air's p-ρ-T data [6] | Not in hand; to be requested and curated in Phase 2 |
| Air density validator | CIPM-2007's compressibility factor of air, Picard et al. [8] | Curated (`picard_2008`): Appendix A1.4, 600–1100 hPa, 15–27 °C, $u(Z) = 1.5 \times 10^{-5}$ |
| Cross pressure virial validator | Bose and Cole, CO2–Ar, $B_{p,12} = -30.1$ cm³/mol at 322.8 K | Curated (`bose_1970`), in notes: to be extracted as a dataset |

**Vocabulary.** Three new quantities, documented on {doc}`../../data` before any data are
added:

- the unlike-pair refractivity virial coefficient (m⁶/mol², abscissa wavelength or
  "static"), with a required parameter `partner`, and the temperature and contribution
  of §7.3;
- the cross second virial coefficient of the pressure (m³/mol, abscissa temperature),
  with the same parameter;
- the third refractivity virial coefficient (m⁹/mol³), for the domain.

A pair is not a material. Each such dataset is filed under the first species of the pair
in sorted order, with the other as `partner`, and the loader makes it symmetric. The
existing `refractivity_virial_coefficient` is not bent to hold a cross term (data rule 7).

**Validators.** The CIPM-2007 air compressibility; Bose and Cole's CO2–Ar mixture pressure
virials; measured cross second virial coefficients for the default pairs. For the
multi-fluid models, the mixture data their papers fit are not independent and are not
used. No independent optical validator of an unlike term is in
hand, and that is stated rather than implied.

**Uncertainty.** A declared unlike term carries the uncertainty of §7.2, and the page
reports its effect on $n - 1$ and on the group-velocity dispersion over the domain. The
mixture density carries its model's stated uncertainty, or, for the fallback, that of
the cross coefficient.

References marked † in §14 are not yet curated source records. Each is checked against
its paper before acceptance, and no value is taken from them until it is.

### 6. Proposed interface

Sketches only, fenced as text until they run.

```text
optmat.gas.mixture.mixture_optics(species) -> MixtureOptics
MixtureOptics.refractivity(wavelength, densities) -> ndarray         # n - 1, real convention
MixtureOptics.susceptibility(wavelength, densities) -> ndarray       # chi, complex
MixtureOptics.dsusceptibility_domega(wavelength, densities) -> ndarray       # s
MixtureOptics.dsusceptibility_ddensity(wavelength, densities) -> dict[str, ndarray]  # m^3
MixtureOptics.pairs -> tuple[PairTerm, ...]      # every (s, t), like and unlike
MixtureOptics.intervals -> tuple[tuple[float, float], ...]   # m, §4.3
MixtureOptics.density_bound(composition) -> float   # m^-3, the optical domain, §4.3
MixtureOptics.identity                           # OEP 1

PairTerm(species, partner, b, sigma, status, source, reason,
         frequency, temperature, contribution)
    # b = B_st / mean(A_s, A_t), m^3/mol, for a curated term; an estimate names its rule
    # status "curated" | "estimate" | "unknown" (a like term optmat has not curated,
    #     and its pairs)
    # frequency "static" or the source's wavelength; temperature in K;
    # contribution "induced" | "total" (a static value that includes orientation)

optmat.gas.mixture.mixture_equation_of_state(composition) -> EquationOfState  # OEP 1's type
EquationOfState.number_density(state), .pressure(number_density, temperature),
    .partial_densities(state)                    # x_s N
EquationOfState.method   # "multi-fluid" | "virial" | "partial-density"
EquationOfState.status   # "qualified" | "provisional"
```

**Real and complex.** `refractivity` keeps the pure models' real-polarizability
convention: for one species it equals `Model.refractivity`. `susceptibility` is complex,
from each constituent's complex polarizability where it has one: for one species it
matches `Model.refractive_index`. The two differ only for an absorbing constituent, by
5 × 10⁻⁸ of $n - 1$ for O3 at 255 nm (R8), and such a constituent warns as its pure model
does.

**Densities.** `densities` maps each species to its number density in m⁻³, so that a
propagation code passes its labelled densities along a fibre directly.
- The wavelength and the densities broadcast together by numpy's rules, and the
  densities may mix scalars and arrays.
- The keys must be exactly the model's species, after aliasing. A missing, extra or
  repeated species is refused by name. Two labelled populations of one species are
  combined by the caller first (§4.1).
- A species at zero density contributes nothing, and all-zero densities give
  $\chi = 0$. The intervals depend on the model's species, not on their densities, so a
  zero density does not change the warnings. This zero-density limit of the optics is
  distinct from OEP 1's refusal of vacuum in the equation of state's inverse.
- The density derivative is taken with respect to number density:
  $\partial\chi/\partial N_s = (1/N_A)\,\partial\chi/\partial\rho_s$, in m³. §4.1 writes
  it in molar density. The frequency dependence enters
through each species' recommended dispersion model (or a model the caller binds);
`wavelength` is the vacuum wavelength in metres, and the range semantics are those of the
constituents' models: the intersection of their qualified intervals (§4.3).
`MixtureOptics.intervals` states it.

For one species, `MixtureOptics` equals the pure model's `refractivity` and
`refractive_index` at the same density; for a composition, `mixture_equation_of_state` equals the pure equation of state
at $x_s = 1$. A species with no dispersion model, or a mixture with a pair the chosen
density method cannot treat, is refused by name.

### 7. Implementation strategies

#### 7.1 Where it belongs

In `optmat.gas`, as `optmat.gas.mixture`: it combines gas models, the gas virials and
the gas equations of state. The core is untouched.

#### 7.2 The unlike-pair refractivity terms

- (a) Zero, with a documented domain (SEP 11's first proposal).
- (b) A declared estimate: the arithmetic mean of the like terms,
  $B_{st} = \tfrac12(B_{ss} + B_{tt})$, each like term dispersing with its own $A$ (§4.1),
  status `estimate`; curated values replace it pair by pair.
- (c) Refuse any mixture without curated unlike terms.

**Decided (2026-09-30): (b), with these uncertainty rules.**
- The standard uncertainty of an estimate is $\max(|B_{ss}|, |B_{tt}|)$. Cancellation
  between like terms of opposite sign, as helium's and argon's, cannot shrink it. It is a
  sensitivity assumption until data give it a statistical basis, and the page says so.
- A like term optmat has not curated (O2 and the others of §2) is `unknown`, not zero.
  It enters with the value zero and a declared scale: $|b| \le 2$ cm³/mol, the largest
  $|b|$ of the curated gases (SF6 and Xe), times the mean $A$. Its pairs are `unknown`
  with the same scale, and their mixtures are provisional.
- Estimated and unknown terms that share a like coefficient are correlated. The mixture's
  sensitivity combines all of them linearly, as fully correlated.

Zero is not a neutral choice: an unlike pair's excess polarizability is of the same kind
and order as a like pair's. The mean of the like terms is defined whatever their signs.
Its uncertainty spans zero and the geometric mean where that exists (for CO2–SF6, 6.57
against a mean of 12.3 cm⁶/mol²). (c) would refuse every mixture, since nothing is
curated.

For the mixtures of §3 the estimate moves $n - 1$ by at most 7 × 10⁻⁵ at 10 bar and
2.2 × 10⁻⁴ at 30 bar. For He/Ar its uncertainty is twice that. It sets the optical bound
of §4.3 for 50/50 He/Ar at about 2700 mol/m³, 65 bar at 294 K.

#### 7.3 The data format

Unlike-pair values are datasets of the new quantities (§5), with `partner`; the choice
per pair, curated or estimated, is recorded in `data/gas/refractivity_virial.toml`
under a `[pairs]` table with its reason, as the like pairs are now. An estimate is a
rule named in that file, not a number typed into it.

Each pair datum states its frequency (static or a wavelength), its temperature and its
contribution: induced, or total including orientation. The selection rules refuse a
static value as an optical term unless its induced part is stated. Bose and Cole's
CO2–Ar $B_{\varepsilon,12}$, whose induced part is assumed, is therefore not an optical
coefficient.

#### 7.4 The mixture density

- (a) The partial-pressure rule: every unlike pressure coefficient set to zero, and the
  composition changed (§3).
- (b) By the evidence for each pair, as below.
- (c) A truncated virial equation for every mixture.

**Decided (2026-09-30): (b).**
- **Multi-fluid** where every pair of the mixture has a binary model fitted to data of
  that pair.
  - Phase 2 compares every method, multi-fluid and fallback, with the measured mixture
    densities and cross second virial coefficients near ambient (§5): 250–350 K, up to
    the optical domain's density. The cross virials were compared but not fitted in
    the models [4]. The densities are data the models fit, and test the fallback
    independently.
  - The default for each mixture is the method closest to the data within their
    uncertainties. The multi-fluid model is `qualified` only where it passes.
  - The acceptance rule is fixed before the comparison: agreement within twice the
    combined uncertainty. A pair that fails is reported, and its mixtures fall back.
  - A placeholder pair (§4.2) counts as unfitted: CoolProp can evaluate it, but it is
    not qualified.
- **The fallback of §4.2** otherwise, with curated cross second virial coefficients for
  every pair without a fitted model: method `virial`, qualified within its truncation
  bound.
- **`partial-density`**, the fallback with every unknown cross coefficient zero, where one
  is missing: status `provisional`, with its stated bound.
- The method and status are part of the equation of state's identity. A provisional one
  warns when it is built.
- **Mixtures required at acceptance:** He/Ar, Ne/Ar, He/Ne and dry air (N2/O2/Ar). None
  of their pairs is a placeholder. Other mixtures run on the fallback, provisional until
  their cross coefficients are curated.

It is the most accurate route where data exist, it reduces to the pure equations exactly,
and it names what is missing where they do not. (c) would discard the pure fluids'
accuracy at the pressures that matter. A propagation code records the method and
status.

For dry air, the multi-fluid N2/O2/Ar mixture and the pseudo-pure fluid OEP 1 uses for
`Gas.Air` differ by 2.0 × 10⁻⁴ at 10 bar. Once this OEP is implemented, the mixture
equation of state serves a labelled air fill and the pseudo-pure fluid serves `Gas.Air`,
and the manual gives their difference.

#### 7.5 Temperature

The like-pair refractivity coefficients are measured at stated temperatures (Achtermann
et al.'s at 323 K). Their temperature dependence is weak and is not modelled: the
coefficients are exposed with their temperatures, as OEP 1 exposes a polarizability's
reference temperature. The pressure coefficients' strong temperature dependence is the
equation of state's, and is modelled.

#### 7.6 Cost and dependencies

The optics are a quadratic form in the densities per frequency: elementwise, as cheap as
the dilute law. A multi-fluid density costs a CoolProp flash, milliseconds per state; a
propagation code evaluates it at the nodes of its density tables, not per step. No new
dependency.

### 8. Verification and validation

| Item | Check | Tolerance set by |
|---|---|---|
| Pure limit | One species equals the pure model's `refractivity` (real) and `refractive_index` (complex) at the same density, for every gas with a model, value and derivative, including an absorbing one (O3 at 255 nm) | Roundoff |
| Interface | Broadcasting, species keys (missing, extra, aliased, repeated), zero and all-zero densities, as §6 states | Exact |
| Invariances | Symmetry in $s, t$; splitting a species into two labels and recombining changes nothing; order of species changes nothing; the dilute limit is $\sum n_s\alpha_s/\varepsilon_0$ | Summation bound |
| Derivatives | Frequency and density derivatives against multi-step central differences | Truncation and rounding, stated |
| Pair terms | Every pair of every mixture has a `PairTerm` with a status; estimates follow the rule of §7.2; `absent` like terms are listed | Exact |
| Sensitivity | The effect of the estimated and unknown terms' uncertainties, combined linearly, on $n - 1$ and GVD over the domain, for He/Ar, Ne/Ar, He/Ne and air, reported on the page | Reported |
| Mixture density | Multi-fluid reduces to each pure fluid at $x_s = 1$. Dry air against CIPM-2007's $Z$ over 600–1100 hPa and 15–27 °C | $u(Z) = 1.5 \times 10^{-5}$, plus the composition difference, computed |
| Composition and inverse | Every method returns $N_s = x_s N$; pressure to density to pressure round-trips, for all three methods | Roundoff and the solver's tolerance |
| Default pairs | Each default pair's effective cross second virial, 250–350 K, against the measured values; every method's density against the measured mixture densities near ambient (§5) | Twice the combined uncertainty, fixed before the comparison; a failure is reported, and the default follows the comparison |
| Cross coefficient | The multi-fluid CO2–Ar model's mixture second virial at Bose and Cole's measured compositions and temperatures (322.8 K; one mixture at 320.8 K) against their mixture pressure virials | Twice the combined uncertainty: their stated 1–3 cm³/mol with the pure coefficients' uncertainties. Fixed before the comparison, reported either way, with Gernert's stated limitation of the model (an effective cross coefficient that depends on composition) |
| Fallback | The fallback against the multi-fluid density on the pairs that have both, within its truncation bound | Measured, reported; checks its truncation assumption |
| Domain | Warnings at the bounds of §4.3, not below them; the third-term bound from the measured $C_R$; `StateError` outside the phase domain | Exact |

### 9. Documentation

- A new manual page, *Gas mixtures at finite density*: the law and its derivation (§4),
  the pair terms and their status per pair, the density methods, the domain, and the
  sensitivity tables; ending with a link to this OEP at Final. Each phase ships its own
  sections of it (§12).
- {doc}`../../conventions`: the mixture law beside the pure one.
- {doc}`../../data`: the two new quantities.
- A gallery example: the dispersion of a He/Ar fill against pressure, separating the
  dilute law, Lorentz–Lorenz, the pair terms and the density rule.

### 10. Compatibility and impact

- **Numbers.** No existing model changes. The pure-gas law is unchanged, and so is every
  fit.
- **Interface.** New names only, in the surface of OEP 1.
- **Register.** Every pair record, estimate rule and validator comparison gets its
  entry.
- **OEP 3.** Its refits of the molecules change the constituents' numbers and intervals.
  The law follows them through the constituents' models, whose digests it binds (OEP 1,
  §4.8).
- **solprop.** SEP 11 consumes the law; the density of mixtures is a change SEP 11 does
  not yet include (its §2 keeps the partial-pressure rule), and §3's numbers suggest it
  should.

### 11. Alternatives considered

- **Per-species Lorentz–Lorenz, then sum.** Wrong at the order of the pair terms (§4.1).
- **Zero unlike terms without an uncertainty** (§7.2(a)). Presents a guess as a value.
- **Keep the partial-pressure density.** Its error is the largest of the finite-density
  effects for these mixtures (§3).
- **A refractometric mixture formula.** None exists for the mixtures of interest beyond
  air, and air's is a different source (OEP 1, §4.4).

### 12. Phases

Each phase is independently mergeable with a green gate, its register entries and its
changelog line.

- **Phase 1 — The optics.**
  - The unlike-pair and third-virial quantities, and the pair record format.
  - Achtermann's $C_R$ extracted.
  - The pair estimates.
  - `MixtureOptics`, with its derivatives and optical domain.
  - The search for unlike-pair data and for the like terms optmat lacks (O2 first), with
    what is found curated.
  - The page's optics sections, and the conventions and data entries.
  - Needs OEP 1's Phase 2b (the derivative).
- **Phase 2 — The density.**
  - The cross second virial quantity; the measured cross virials and mixture densities
    for the default pairs, and the cross virials of the fallback pairs solprop uses;
    Bose and Cole's mixture virials extracted.
  - `mixture_equation_of_state` with its three methods, statuses, phase checks and
    validators; the CoolProp pair records.
  - The page's density sections.
  - Needs OEP 1's Phase 2c (the equation-of-state object).
- **Phase 3 — The example and the combined sensitivity discussion.**

The roadmap's informal item on the finite-density mixture law (its *solprop* section) is
covered in full.

### 13. Open questions

- **Q1.** The estimate rule of §7.2: the arithmetic mean with a 100 % uncertainty, or a
  rule built on the pair polarizability's long-range form? The recommendation is the
  simple rule until data arrive.
  *Decided (2026-09-30):* the mean of the like terms' $B$, with the uncertainty rules of
  §7.2.
- **Q2.** Which mixtures does solprop use, beyond He/Ar, Ne/Ar and air? The pairs without
  multi-fluid parameters (He–Xe, Ne–Xe, Ar–Xe and the SF6 and N2O pairs) need cross
  second virial coefficients for the fallback, and those are to be requested.
  *Decided (2026-09-30):* He/Ar, Ne/Ar, He/Ne and dry air at acceptance (§7.4). Other
  mixtures run on the fallback, provisional until their cross coefficients are curated.
- **Q3.** The domain's 1 % criterion (§4.3): is 10⁻⁴ of $y$ the right bound for the
  neglected third term?
  *Decided (2026-09-30):* 10⁻³ of $y$, bounded with the measured $C_R$ (§4.3).
- **Q4.** Should O2's like-pair refractivity coefficient be sought before Phase 1, since
  air's optics depend on it and on N2–O2?
  *Decided (2026-09-30):* in Phase 1's search, not before acceptance. Until then O2 is
  `unknown` under §7.2, and air's sensitivity to it is reported.

The following questions record the GPT-6 review of 2026-09-30. Each now carries its
disposition (2026-09-30).

- **Q5.** What evidence bounds the omitted third refractivity term? The 1 % rule does
  not establish the claimed error; see {ref}`R1 <oep4-review-r1>` and Q3.
  *Resolved in the text:* Achtermann's measured $C_R$, a declared provisional scale where
  none is measured, and no squared-ratio argument (§4.3).
- **Q6.** How are missing coefficients, cancellation and uncertainty distinguished from
  measured zero coefficients? See {ref}`R2 <oep4-review-r2>`, Q1 and Q4.
  *Resolved in the text:* the status `unknown` with a declared scale, an uncertainty that
  cancellation cannot shrink, and a linear combination of shared terms (§7.2).
- **Q7.** What is the complete fallback equation, how does it preserve the declared
  composition, and how is its inverse verified? See {ref}`R3 <oep4-review-r3>`.
  *Decided (2026-09-30):* the partial-density law with the cross second virial (§4.2),
  which preserves the composition exactly. The round trip is tested for every method
  (§8).
- **Q8.** Does the estimate average $B$ or $b = B/A$, and which regenerated numerical
  evidence supports that choice? See {ref}`R4 <oep4-review-r4>` and Q1.
  *Decided (2026-09-30):* $B$, which §3's table and the CO2–SF6 example already use. §3
  now states the references, the air composition and that the columns are relative
  changes.
- **Q9.** Which temperature, wavelength and equation-of-state inputs determine the
  optical density bound and the mixture phase domain? See {ref}`R5 <oep4-review-r5>`.
  *Resolved in the text:* the optical bound needs no temperature and holds across the
  intervals; phase admissibility belongs to the equation of state, which raises
  `StateError` (§4.3).
- **Q10.** What qualifies a supported backend mixture as a recommended model, and over
  what composition and state range? See {ref}`R6 <oep4-review-r6>`.
  *Decided (2026-09-30):* a binary model fitted to data of the pair, and validated near
  ambient against measured mixture densities and cross second virials, in comparison
  with the fallback. Placeholder pairs are not qualified (§7.4).
- **Q11.** What observable and uncertainty budget define the CO2–Ar cross-virial
  acceptance test? See {ref}`R7 <oep4-review-r7>`.
  *Resolved in the text:* the mixture second virial at Bose and Cole's compositions, with
  the pure coefficients' uncertainties, under a rule fixed before the comparison (§8).
- **Q12.** Which real or complex refractivity convention and absorbing-gas pure limit
  does the interface promise? See {ref}`R8 <oep4-review-r8>`.
  *Decided (2026-09-30):* both. `refractivity` is real, with the pure limit
  `Model.refractivity`, and `susceptibility` is complex, with the pure limit
  `Model.refractive_index`. O3 is the absorbing test case (§6, §8).
- **Q13.** How is the molar-density derivative converted to the public number-density
  derivative? See {ref}`C1 <oep4-review-comments>`.
  *Resolved in the text* (§6).
- **Q14.** What are the broadcasting, species-key and zero-density rules, including
  their effect on range warnings? See {ref}`C2 <oep4-review-comments>`.
  *Resolved in the text* (§6, §8).
- **Q15.** Which required metadata prevent a static molecular dielectric cross term
  from being used as an optical coefficient? See {ref}`C3 <oep4-review-comments>`.
  *Resolved in the text:* frequency, temperature and contribution on every pair datum
  (§6, §7.3).
- **Q16.** Which manual documentation ships with each independently mergeable phase?
  See {ref}`C4 <oep4-review-comments>`.
  *Resolved in the text* (§9, §12).
- **Q17.** Do the declared scales `unknown_b` and `third_ratio` change, given xenon's b of
  2.46 and hydrogen's unresolved ratio of 1679 (TC-774)? Raised in Phase 1.
  *Decided (2026-10-01):* `unknown_b` becomes 2.5 cm³/mol, the largest curated |b| as
  §7.2 words it; `third_ratio` stays 600, with H2 excluded because its B_R is not
  resolved (TC-795).
- **Q18.** Does neon's B_R (with krypton's and xenon's) move to Achtermann et al. 1993,
  and does oxygen take Egan et al.'s b_R (TC-771, TC-775)? Raised in Phase 1.
  *Decided (2026-10-01):* neon and krypton take Achtermann et al. 1993's optical values
  (TC-792, TC-793); xenon keeps Burns et al.'s; oxygen takes Egan et al.'s estimate as a
  provisional like term with its 50 % uncertainty (TC-794).
- **Q19.** Are the density defaults of Phase 2 accepted, and should air's range stay at
  Martin et al.'s 290–320 K (TC-779)? Raised in Phase 2.
  *Decided (2026-10-01):* the defaults of TC-779 are accepted for the noble-gas
  mixtures; dry air's default is the multi-fluid model at every temperature,
  provisional, over the rule of TC-778 (TC-796).

### 14. References

1. A. D. Buckingham and J. A. Pople, "The dielectric constant of an imperfect non-polar
   gas," *Trans. Faraday Soc.* **51**, 1029 (1955). doi:10.1039/tf9555101029
   (`buckingham_1955`)
2. † J. H. Dymond, K. N. Marsh and R. C. Wilhoit, *Virial Coefficients of Pure Gases and
   Mixtures: Virial Coefficients of Mixtures*, Landolt–Börnstein, Group IV (Springer),
   introduction. doi:10.1007/10754889_1. Could not be obtained (2026-09-30); nothing is
   taken from it.
3. O. Kunz and W. Wagner, "The GERG-2008 wide-range equation of state for natural gases
   and other mixtures: an expansion of GERG-2004," *J. Chem. Eng. Data* **57**, 3032–3091
   (2012). doi:10.1021/je300655b (`kunz_2012`)
4. J. Tkaczuk, I. H. Bell, E. W. Lemmon, N. Luchier and F. Millet, "Equations of state
   for the thermodynamic properties of binary mixtures for helium-4, neon, and argon,"
   *J. Phys. Chem. Ref. Data* **49**, 023101 (2020). doi:10.1063/1.5142275
   (`tkaczuk_2020`)
5. G. Garberoglio and A. H. Harvey, "Path-integral calculation of the second dielectric
   and refractivity virial coefficients of helium, neon, and argon," *J. Res. NIST*
   **125**, 125022 (2020). doi:10.6028/jres.125.022 (`garberoglio_2020`; cited by
   solprop's SEP 11)
6. E. W. Lemmon, R. T. Jacobsen, S. G. Penoncello and D. G. Friend, "Thermodynamic
   properties of air and mixtures of nitrogen, argon, and oxygen from 60 to 2000 K at
   pressures to 2000 MPa," *J. Phys. Chem. Ref. Data* **29**, 331–385 (2000).
   doi:10.1063/1.1285884 (`lemmon_2000`)
7. † J. H. Dymond, K. N. Marsh and R. C. Wilhoit, as [2], the tables of mixtures
   (chapters doi:10.1007/10754889_n). Could not be obtained; the primary papers of §5
   replace it.
8. A. Picard, R. S. Davis, M. Gläser and K. Fujii, "Revised formula for the density of
   moist air (CIPM-2007)," *Metrologia* **45**, 149 (2008).
   doi:10.1088/0026-1394/45/2/004
9. G. J. Gernert, *A New Helmholtz Energy Model for Humid Gases and CCS Mixtures*
   (doctoral thesis, Ruhr-Universität Bochum, 2013). A public copy is hosted by NIST
   (`gernert_2013`).

## Part II — Implementation record

### Phase 1 — The optics

- **Scope:** Phase 1 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0004-optics`, merged into `main`. No recommended value changes.
- **Delivered:**
  - `optmat.gas.mixture`:
    - `mixture_optics(species)` and `MixtureOptics`, with `refractivity` (real),
      `susceptibility` (complex), `dsusceptibility_domega`, `dsusceptibility_ddensity`,
      `refractivity_uncertainty`, `pairs`, `intervals`, `density_bound(composition)`
      and the identity;
    - `PairTerm`.
  - `optmat.gas.virial`: `like_term` (with σ, frequency, temperature and contribution),
    `third_term` and `pair_rules`.
  - `data/gas/refractivity_virial.toml`: the `[pairs]` rules, the `C_R` pointers, helium's
    `C_R_bound`, and the contribution of the static values.
  - The quantity `unlike_pair_refractivity_virial_coefficient` (no datasets).
  - The data:
    - Achtermann et al. 1991's C_R and D_R, extracted;
    - Achtermann et al. 1993 and Huot and Bose 1991, curated.
  - The page {doc}`../../mixtures`, with the sensitivity table of §8; the mixture law on
    {doc}`../../conventions`; the quantity on {doc}`../../data`.
  - `tests/test_mixture.py`; new loader tests in `tests/test_virial.py`; 8 golden
    records.
- **Deviations from the proposal:**
  - C_R needs no new quantity: `refractivity_virial_coefficient` with `order = 3` holds
    it (§5 lists it as new) (TC-772).
  - A repeated species key is not refused, because it cannot arise: a `Gas` and its name
    are one dictionary key (TC-773).
  - The density bound is the smallest over 48 wavelengths per qualified interval and the
    static limit. The figure in §7.2 is at one wavelength: 50/50 He/Ar is 2446 mol/m³
    here, against about 2700.
  - §6's "5 × 10⁻⁸ of $n - 1$" for O3 is an absolute difference, −5.27 × 10⁻⁸, about
    1.2 × 10⁻⁴ of $n - 1$; it is $-(\operatorname{Im}\chi)^2/24$ at second order (TC-773).
  - A mixture with only curated and estimated terms is `validated`, as §7.2's wording
    implies. No unlike term has been checked against a measurement, and the page says so.
- **Decisions made:**
  - A static B_R is an optical term only where `contribution = "induced"` is stated; Ne
    and Kr, atoms, state it (TC-772).
  - A gas with neither B_R nor C_R gives no third-term bound; its unknown like term
    bounds it instead, at 500 mol/m³ alone (TC-773).
  - The uncertainty of an unknown pair is the larger of its like terms', the unknown one
    at its scale. This is §7.2's estimate rule applied to an unknown like term.
- **Register entries:** TC-771 to TC-775.
- **Verification:**
  - The full gate, `uv run poe check`, on the branch before the merge.
  - The library test for the three extractions (`achtermann_1991`, `achtermann_1993`,
    `huot_1991`).
  - The pure limit for all 15 gases:
    - value within 1e-13 of `Model.refractivity`;
    - complex index within 1e-15 of `Model.refractive_index`;
    - frequency derivative within 1e-13 of the pure law's chain rule.
  - O3 at 255 nm reproduces R8's four printed figures.
  - Both derivatives agree with fourth-order differences.
  - The interface rules of §6 and the invariances of §8.
  - The unlike term and the sensitivity are checked against direct computation.
  - The third-term bounds of argon and N2 against §4.3.
  - Not run: `test-fits`, since no fit changed.
- **Sources checked:**
  - Against the page images: Achtermann et al. 1993's Table II, Huot and Bose's Table I, and
    Achtermann et al. 1991's C_R and D_R columns.
  - Against these values, neon's comparison of TC-744 is confirmed (TC-771).
  - The search for unlike-pair data and O2's B_R found none for the pairs of interest.
    Egan et al.'s O2 b_R, a heuristic with 50 % uncertainty, is the one candidate
    (TC-775).
- **Limitations:**
  - Every unlike term is an estimate, and no optical validator of one is in hand.
  - The like terms' dispersion is the pure law's convention, and their temperature
    dependence is neglected (§7.5).
  - Air's domain, 1296 mol/m³, is set by oxygen's unknown like term.
- **Open questions:**
  - Q17. The declared scales against the curated coefficients (TC-774). Xe's b, 2.46,
    exceeds `unknown_b = 2` which §7.2 calls the largest curated; H2's |C_R| A_R/B_R², 1679,
    exceeds §4.3's "largest measured ratio, N2's 560", though its B_R is not resolved.
    Should either scale change? For the author.
  - Q18. Neon's B_R (TC-744, TC-771). Should it, and possibly krypton's and xenon's, move
    to Achtermann et al. 1993? Should oxygen take Egan et al.'s b_R (TC-775)? These change
    recommended values, so they are prepared on a separate branch for the author.
- **Handoff:**
  - Phase 2 adds `mixture_equation_of_state`, which gives the number densities that
    `MixtureOptics` takes.
  - Curated unlike-pair data, when found, are datasets of the new quantity, chosen per
    pair under `[pairs]`, and replace the estimate pair by pair.

### Phase 2 — The density

- **Scope:** Phase 2 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0004-density`, merged into `main`. New functionality only; no
  existing value changes.
- **Delivered:**
  - `optmat.gas.mixture_eos`:
    - `mixture_equation_of_state(composition, method=None)` and
      `MixtureEquationOfState`, with `number_density`, `partial_densities`, `pressure`,
      `compressibility`, `second_virial`, `pairs`, `method`, `qualified`, the phase
      checks and the identity;
    - `PairDensity`.
  - `data/gas/mixture_pairs.toml`, the curated pairs. Each model key is checked against
    CoolProp's, with its qualification, the fallback's cross coefficient, and the
    fallback's qualification.
  - `optmat.gas.mixture_validation`: the comparison of §8 under the rule fixed
    beforehand (TC-778).
  - The quantities: a `thermodynamic` family of six (pure, cross, excess and mixture
    second virials; mixture compressibility factor and density), and new pressure and
    molar-density units.
  - The data: seven new records (Brewer and Vaughn, Martin et al., Vatter et al., Vogl
    and Hall, Kalfoglou and Miller, Blancett et al., Kestin et al.), and Bose and
    Cole's pressure virials.
  - The density sections of {doc}`../../mixtures`, with its snippets now run by the docs
    test; `tests/test_mixture_eos.py` and `tests/test_mixture_validation.py`; 5 golden
    records.
- **Deviations from the proposal:**
  - The function lives in `optmat.gas.mixture_eos`, not `optmat.gas.mixture`. It returns
    `MixtureEquationOfState`, which has OEP 1's method names but is not its
    `EquationOfState`.
  - Its status is the package's `validated` or `provisional`, with a `qualified`
    field, in place of a "qualified" status (TC-785).
  - The fallback's cross coefficient may be a measured excess $E$, with
    $B_{12} = E + (B_{11} + B_{22})/2$ from the reference equations' pure $B$. This was
    amended after the comparison, as TC-779 records with the reason.
  - The same pass rule is applied to the fallback against the measured densities. It
    can only make a default provisional.
  - The fallback holds only at its cross coefficients' measured temperatures.
  - The CO2–Ar check of §8 is reported (TC-781), not used for qualification: Bose and
    Cole's $B_{p,12}$ states no uncertainty.
- **Decisions made:**
  - The acceptance rule and the default order, fixed before the comparison (TC-778,
    amended before it to count $E$).
  - The pair records follow the comparison, and a test keeps them in step.
- **Results** (TC-779):
  - No required pair's multi-fluid model qualifies. Tkaczuk et al.'s noble-gas models
    miss the measured cross virials by 3 to 11 cm³/mol and the densities by up to
    1 × 10⁻² in Z. CoolProp reproduces their own check values, so the gap is the
    models'.
  - Gernert's air pairs miss Martin et al.'s excess coefficients by about
    0.5 cm³/mol.
  - The defaults:
    - He/Ne: the fallback, validated (35 of 36 of Vogl and Hall's points);
    - He/Ar, Ne/Ar and dry air: the fallback, provisional;
    - dry air at 290–320 K: within 1.8 × 10⁻⁶ of CIPM-2007's $Z$.
  - Ar–CH4, CH4–CO2 and CH4–N2 (Kunz and Wagner) qualify.
- **Register entries:** TC-778 to TC-785.
- **Verification:**
  - The full gate, `uv run poe check`, on the branch before the merge.
  - The library test for the eight extractions.
  - Every method preserves the composition and round-trips within 1e-12.
  - The pure limits hold within 1e-13 (multi-fluid) and 1e-10 (fallback).
  - The fallback is the pressure sum of §4.2.
  - The ratio to the partial-density rule follows §4.2's first-order law.
  - CIPM-2007's dry air, and the pseudo-pure air's 2.0 × 10⁻⁴ of §7.4.
  - Tkaczuk et al.'s Table 8.
  - The phase refusals, and the provisional and truncation warnings.
- **Sources checked:**
  - Each new record's tables against the page images, with the reproductions their
    register entries give.
  - The mixture models' parameters against Tkaczuk et al.'s Table 6, and their check
    values against Table 8.
  - CIPM-2007's $Z$ against picard_2008's recorded formula.
- **Limitations:**
  - The noble-gas defaults are provisional. The fallback leaves out the unlike third
    virials and misses Blancett et al.'s helium–argon data by up to 7 × 10⁻⁴ in $Z$ at
    28 bar.
  - Air's default holds only at 290–320 K, and its pairs have no density data of
    their own.
  - Blancett et al.'s scan is incomplete (TC-784).
  - Kalfoglou and Miller's mixture data need a 1/atm unit (TC-783).
  - The mixture phase check of the fallback is a necessary condition only.
- **Open questions:**
  - Q19. Are these defaults accepted? Should air's range stay at 290–320 K or take the
    multi-fluid model outside it? Are more density data sought for Ne–Ar and air's
    pairs, and a complete Blancett et al.? For the author (TC-779).
- **Handoff:**
  - Phase 3 is the He/Ar gallery example and the combined sensitivity tables. It uses
    `partial_densities` into `MixtureOptics`.

### Phase 3 — The example and the combined sensitivity

- **Scope:** Phase 3 of §12, 2026-10-01, implemented by Opus 5.5 (AI) for John Travers,
  on `feature/oep-0004-example`, merged into `main`.
- **Delivered:**
  - the gallery example *The dispersion of a helium-argon fill against pressure*:
    the dilute law, the pair terms, the density rule and the density method, each
    against the full law, in $n - 1$ and in the GVD, with the unlike estimate's band;
  - the section *How much each piece matters* on {doc}`../../mixtures`, whose table a
    test recomputes.
- **Deviations from the proposal:** the example adds the density method's spread
  (multi-fluid against the default fallback) to §9's four pieces, since Phase 2 found
  it the second-largest.
- **Results (TC-786):**
  - At 10 bar the partial-pressure density changes $n - 1$ by 6 × 10⁻³ and the density
    method by 2 × 10⁻³.
  - The local field changes it by 1.1 × 10⁻³, and the pair terms and their declared
    uncertainty by 1.5 × 10⁻⁴ each.
  - The density limits a helium-argon fill.
- **Register entries:** TC-786.
- **Verification:**
  - The full gate, with the gallery built.
  - The page's table recomputed by `tests/test_mixture.py`, both $n - 1$ and the GVD,
    to the printed figures.
- **Sources checked:** none new; the example uses the curated data of Phases 1 and 2.
- **Limitations:** one mixture, one wavelength and one temperature. The tables of
  Phase 1 (He/Ar, Ne/Ar, He/Ne and air) give the other mixtures' optical
  sensitivities.
- **Open questions:** none new.
  - Q17-Q19 remain for the author: the declared scales, the like-pair B_R choices, and
    the density defaults.
  - Every phase is merged and reported. The Final status is the author's to approve
    (OEP 0), once Q17-Q19 are settled or moved to the roadmap.
- **Handoff:** solprop's SEP 11 can take `mixture_equation_of_state(...)
  .partial_densities(state)` into `mixture_optics(...)`. It should record the
  equation's method and status, and its provisional warnings.

## Part III — History

- 2026-09-30: Drafted, at solprop SEP 11's request for the upstream mixture law. The
  sizes of §3 and §4.1 were computed with optmat at `ea0dcaf` and CoolProp 8.0.0; the
  unlike-pair data curated in optmat were read from their source records; the
  bibliographic data of the † references were checked against their DOI records, and
  the papers themselves were not read.
- 2026-09-30: The wavelength domain of a mixture added (§4.3, §6, §10): the intersection
  of its constituents' qualified intervals (OEP 1), with a rule weighted by mole fraction
  to come from OEP 3.
- 2026-09-30: GPT-6 (AI) reviewed the whole of Part I against OEP 0, the data
  conventions, OEPs 1 and 3, solprop SEP 11 and the current implementations. Changes
  requested before acceptance: eight substantive findings and four further comments,
  recorded below with their evidence and as open questions Q5–Q16. The proposal remains
  Draft; no human review or acceptance is implied by this record.
- 2026-09-30: Aligned with OEP 3's revision for its review (§4.3). A mixture widens the
  intersection of its constituents' intervals only under OEP 3's three checks. A
  constituent's different representation, needed to cross its kept bands, is a
  different model.
- 2026-09-30: John Travers (human) reviewed Part I. With GPT-6's review, the review
  status is `Reviewed`. Q1-Q16 remain open before acceptance.
- 2026-09-30: Resolved Q1–Q16 (Opus 5.5 drafting; the author's decisions on Q1–Q4, Q7,
  Q8, Q10 and Q12).
  - The estimate is the mean of the like terms' $B$. Its uncertainty cannot cancel, and
    an absent like term is `unknown` with a declared scale.
  - The optical domain bounds the third term at 10⁻³ of $y$ with Achtermann's measured
    $C_R$, and the declared-term uncertainty at 10⁻³ of $n - 1$. Phase admissibility
    belongs to the equation of state.
  - The fallback is the partial-density law with the cross second virial, which
    preserves the composition.
  - Multi-fluid is the default only for pairs fitted to their own data and validated
    against measured cross second virials. He/Ar, Ne/Ar, He/Ne and air are required at
    acceptance.
  - `refractivity` keeps the real convention, and `susceptibility` is complex.
  - The interface's units, broadcasting, keys and zero densities are stated, pair data
    carry frequency, temperature and contribution, and each phase ships its manual
    sections.
  - §3 now states its references and composition, and adds the partial-pressure rule's
    change of composition.
  - The review record below is kept unchanged as evidence.
- 2026-09-30: Checked against the curated papers (TC-738 to TC-743).
  - §4.1's pair law is Buckingham and Pople's, factors included.
  - Quantum corrections for He and Ne are below 1 % at room temperature.
  - CoolProp's air pairs are GERG-2008's and Gernert's, not Lemmon et al.'s mixture
    model.
  - Near ambient, Tkaczuk et al.'s models deviate from measured He–Ar and Ne–Ar
    densities by more than their difference from the partial-pressure rule. Phase 2
    therefore compares every density method with measured densities and cross virials.
  - The quadratic virial law is cited to Gernert, since the Landolt–Börnstein
    compilation could not be obtained; primary papers replace it for the validation
    data.
  - Gernert's name and Tkaczuk et al.'s article number are corrected.
- 2026-09-30: Accepted by John Travers. Part I is frozen.

(oep4-review-20260930)=
### Review of 2026-09-30 — GPT-6 (AI)

**Recommendation: changes requested before acceptance.** The central structure is
sound: sum the species and pair contributions to the Lorentz–Lorenz function and
invert once. The accuracy claims, numerical evidence and public contract need the
changes identified below. P1 identifies a finding that undermines a central accuracy
or composition guarantee; P2 identifies a numerical inconsistency or an incomplete
contract or validation rule that also needs resolution before acceptance.

**Scope and baseline.** The review read the complete proposal at optmat main commit
`ba368030570ff3f88322f6a7aabfb3c6dc9911df`, alongside OEPs 0, 1 and 3, solprop's SEP 11
and `solprop.physdata.mixtures`, optmat's conventions, source records, virial choices,
equation of state and gas model implementations. Numerical probes used CoolProp 8.0.0
and the baseline's recommended optical models. They are review calculations, not new
curated datasets or changes to the package's models. References to proposal sections
below concern that reviewed version; later resolutions should state what changed.

(oep4-review-r1)=
#### R1 (P1) — The density cutoff does not bound the omitted third term

**Where:** §4.3's first density criterion; §4.1's neglected third term; Q3.

Requiring the second contribution to remain below 1 % of the linear contribution
does not establish that the third is below $10^{-4}$. Write

$$
r_2 = \frac{B_R\rho}{A_R}, \qquad
r_3 = \frac{C_R\rho^2}{A_R}
= \frac{C_R A_R}{B_R^2}\,r_2^2.
$$

The squared-ratio argument requires a bound on $C_R A_R/B_R^2$. Achtermann et al.,
Table II, p. 5675 (`achtermann_1991`, doi:10.1063/1.460478), give nitrogen's values
at 323 K and 632.99 nm as $A_R = 4.4464$ cm³/mol, $B_R = 0.89$ cm⁶/mol² and
$C_R = -99.8$ cm⁹/mol³. Their ratio is approximately −560.22, not of order unity.
The source record already records $C_R$ in its notes. The reviewer checked the
original table image, including its units and temperature.

Using `optmat.gas.eos.number_density` at 323 K and 100 bar gives the following
counterexamples, evaluated with the printed optical coefficients:

| Gas | Molar density (mol/m³) | $B_R\rho/A_R$ | $C_R\rho^2/A_R$ |
|---|---:|---:|---:|
| N2 | 3660.813 | $7.3276\times10^{-4}$ | $-3.0080\times10^{-4}$ |
| Ar | 3825.732 | $1.6047\times10^{-3}$ | $-3.0207\times10^{-4}$ |

Argon's coefficients on the same page are $A_R = 4.196$ cm³/mol,
$B_R = 1.76$ cm⁶/mol² and $C_R = -86.6$ cm⁹/mol³. Both examples pass the 1 %
second-term criterion but exceed the claimed third-term bound. They are
supercritical states, so excluding condensed pure phases does not remove the
counterexamples. In a mixture, cancellation between positive and negative pair
contributions can remove the second-term cutoff without removing the third term.

**Requested resolution:** use measured higher coefficients to establish bounds where
available; elsewhere declare a provisional domain with an explicit assumption or
validation. Do not present the square of the second-term ratio as an error bound.
This can be done while leaving the third term out of the evaluator.
**Status:** open, Q5.

(oep4-review-r2)=
#### R2 (P1) — Missing coefficients can produce zero estimated uncertainty

**Where:** §7.2's estimate and absent-term rules; §4.3's uncertainty cutoff; Q1 and Q4.

The proposed uncertainty equals the magnitude of the estimate, while an unavailable
like coefficient enters as zero. Applied literally to O2–D2, whose like terms are
both absent, it gives $b_{st} = 0$ and an uncertainty of zero. Opposite-sign like
coefficients can also make the mean, and therefore its assigned uncertainty, small
through cancellation. Neither case establishes that the unlike interaction is small.

An `absent` status records provenance but does not supply the uncertainty needed to
compute `density_bound`. If inheriting that status instead means an unavailable
uncertainty, the proposal needs to say what the bound does in that case. Air's reported
sensitivity also needs to account for its missing O2 self term, not only the estimated
unlike terms.

**Requested resolution:** distinguish unknown from measured zero, and supply an
independently justified uncertainty scale or explicitly leave the uncertainty and
resulting qualification unavailable. Until it has a statistical justification, the
±100 % range is a sensitivity assumption rather than a standard uncertainty. Record
how uncertainties shared through the same like coefficients are combined.
**Status:** open, Q6.

(oep4-review-r3)=
#### R3 (P1) — The fallback does not define composition preservation or a full inverse

**Where:** §4.2, §6's `partial_densities` and inverse, and §7.4's fallback.

The partial-pressure rule evaluates

$$
N_s = \frac{x_s p}{Z_s(x_s p,T)\,k_B T}.
$$

Its resulting number fractions $N_s/\sum_t N_t$ generally differ from the input mole
fractions $x_s$. A correction to the total density alone does not resolve this.
For nominal 50/50 He–Ar at 294 K, the review computed:

| Pressure | He fraction from the partial-pressure evaluations | Relative change of $n-1$ against the same optical law at the multi-fluid density and fixed 50/50 composition |
|---|---:|---:|
| 10 bar | 0.49856236 | $+3.8472\times10^{-3}$ |
| 30 bar | 0.49573890 | $+1.2208\times10^{-2}$ |

The optical comparison uses 800 nm, the same like coefficients and the proposed
mean-$b$ unlike rule in both evaluations; only the constituent densities change.
The differences therefore exceed §3's total-density differences. They are comparisons
of two model constructions, not independently established experimental errors.

The fallback is also given only as a low-density expansion. For example, with
$\Delta B = \sum_{s\ne t}x_sx_t B_{st}^{(p)}$, both
$\rho = \rho_{\rm pp}(1-\Delta B\rho_{\rm pp})$ and
$\rho = \rho_{\rm pp}/(1+\Delta B\rho_{\rm pp})$ reproduce the stated leading
correction. They differ at higher orders and have different inverses. Neither is
selected by the proposal. The final provisional fallback also needs a rule for its
stated bound when the cross coefficients are unavailable.

**Requested resolution:** define one complete forward equation and its inverse, with
$N_s = x_s N_{\rm total}$ for an EOS bound to a fixed composition. Specify the
temperature and density range of the correction, including its truncation uncertainty,
and add composition-preservation and pressure–density round-trip checks for all three
methods. Reconcile the changed interpretation with solprop's labelled fills.
**Status:** open, Q7.

(oep4-review-r4)=
#### R4 (P2) — The size table and recommendation average different quantities

**Where:** §3's size table, §4.1's dispersion convention and §7.2's estimate.

The recommended rule averages $b = B/A$, so that

$$
B_{st}^{\text{mean }b}
= \frac{b_{ss}+b_{tt}}{2}\,\frac{A_s+A_t}{2}.
$$

The numbers in §3 instead match $B_{st}^{\text{mean }B} = (B_{ss}+B_{tt})/2$.
The distinction is substantial for species with different polarizabilities. At
800 nm, 294 K and 10 bar, adding the unlike term gives:

| Mixture | Relative change of $n-1$, recommended mean-$b$ rule | Relative change, mean-$B$ rule (approximately) | §3 |
|---|---:|---:|---:|
| He/Ar 50/50 | $2.95\times10^{-5}$ | $7.33\times10^{-5}$ | $7\times10^{-5}$ |
| Ne/Ar 50/50 | $1.22\times10^{-5}$ | $5.74\times10^{-5}$ | $6\times10^{-5}$ |

The reference is the same mixture at the same multi-fluid density, retaining its like
terms but setting the unlike term to zero. The mean-$B$ figures above use the
leading-order relative contribution; the recommended mean-$b$ figures use the full
inversion. At 30 bar the full mean-$b$ changes are $8.83\times10^{-5}$ for He/Ar
and $3.65\times10^{-5}$ for Ne/Ar, also inconsistent with the table's
$2.2\times10^{-4}$ and $1.7\times10^{-4}$.

The CO2–SF6 example's 12.3 cm⁶/mol² is likewise the arithmetic mean of the printed
$B_R$ values 1.9 and 22.7. Applying the proposed rule with the cited molar
refractivities 6.644 and 11.323 cm³/mol instead gives approximately
10.29 cm⁶/mol². These inputs are recorded in `hourri_1997`'s notes as reproduced
from Achtermann. The review did not freshly verify Hourri's original page.

**Requested resolution:** choose and justify the quantity being averaged, then
regenerate the table and uncertainty comparisons. State the air composition and
normalisation, the reference in each comparison and that the optical columns are
relative changes, not absolute increments of $n-1$.
**Status:** open, Q8.

(oep4-review-r5)=
#### R5 (P2) — The domain interface lacks required inputs

**Where:** §4.3's domain, §6's `density_bound(composition)` and the phase dependencies
in §12.

Phase limits depend on temperature and the chosen EOS. For the proposed cross-term
dispersion, the mixture's pair-to-linear ratio generally depends on wavelength as
well. `density_bound(composition)` takes none of these inputs, and
`mixture_optics(species)` binds neither a temperature nor an EOS.

The minimum of the pure-component phase domains is also not a definition of the
mixture's phase stability. The mixture has a composition-dependent phase envelope.
A pure-fluid saturation limit can reject a permissible dilute constituent in a gas
mixture; the pure limits do not establish that a mixture cannot separate into phases.
CoolProp's documented mixture PT flash performs mixture stability and phase-split
calculations, rather than inferring stability from the constituents' pure phases.

**Requested resolution:** separate optical qualification from thermodynamic phase
admissibility, and specify the inputs or construction-time bindings of each. State
whether the optical bound applies at a wavelength or conservatively across stated
intervals. Specify when invalid phase states raise `StateError` rather than an optical
`ModelRangeWarning`, consistently with OEP 1. Reconcile Phase 1's promised domain with
Phase 2 delivering the mixture EOS.
**Status:** open, Q9.

(oep4-review-r6)=
#### R6 (P2) — Backend support does not establish model qualification

**Where:** §3's density-error claims, §4.2's supported pairs and §7.4's recommended
selection rule.

The review reproduced 39 supported pairs among the fourteen gases. Their backend
reference keys comprise 17 `Kunz-JCED-2012`, ten `Gernert-Thesis-2013`, nine
`Bell-JCED-2016` and three `Tkaczuk-JPCRD-2020` records. Parameter availability
therefore covers different models and levels of evidence. He–H2, He–O2, He–CO,
Ar–H2, H2–O2 and O2–CO all have the four reducing parameters equal to one and
$F_{ij} = 0$, with no departure contribution. Existence of such a record does not
establish that measurements of that pair qualified its prediction.

Tkaczuk et al. [4], §3.1, report percent-level density deviations over their
He–Ar and Ne–Ar datasets. Their discussion gives 0.5–2.5 % for 95 % of the He–Ar
data in the 0–10 MPa comparison and 3 % for 90 % of the Ne–Ar density data. These
figures are not pointwise uncertainties at 294 K, and should not be applied as such.
They nevertheless prevent a backend comparison alone from establishing a
$10^{-4}$-level physical error of the partial-pressure rule.

**Requested resolution:** qualify each recommended EOS by the evidence for its pair
parameters, composition and state range, and uncertainty. Define how unsupported
qualification differs from unsupported evaluation. Describe §3's results as differences
between models until independent evidence establishes which is closer to experiment;
carry the model uncertainty into the comparison and default-selection rationale.
**Status:** open, Q10.

(oep4-review-r7)=
#### R7 (P2) — The CO2–Ar validator needs an observable and an uncertainty budget

**Where:** §5's cross pressure validator and §8's cross-coefficient acceptance check.

For a binary mixture with CO2 fraction $x$, the reviewer computed the effective cross
coefficient from the backend's low-density virial as

$$
B_{12}^{\rm eff}(x,T)
= \frac{B_{\rm mix}(x,T)-x^2B_{11}(T)-(1-x)^2B_{22}(T)}{2x(1-x)}.
$$

CoolProp 8.0.0's `HEOS` CO2–Ar model uses `Gernert-Thesis-2013`. Evaluating its
`Bvirial` and the pure-fluid virials at 322.8 K gives:

| CO2 mole fraction | Effective $B_{12}$ (cm³/mol) |
|---|---:|
| 0.1 | −36.52 |
| 0.5 | −34.11 |
| 0.9 | −34.56 |

These are virial-limit coefficients; the probe initialised each state at 1000 Pa and
used the backend's `Bvirial` evaluator, not a finite-pressure estimate of $(Z-1)/\rho$.
The model therefore does not supply a unique composition-independent coefficient
matching the exact quadratic law. Gernert discusses this limitation of the
corresponding-states mixture formulation in §5.2, pp. 128–130 of *A New Helmholtz
Energy Model for Humid Gases and CCS Mixtures* (Ruhr-Universität Bochum, 2013;
[public NIST copy](https://trc.nist.gov/refprop/Documents/Humid_Gases_and_CCS_Mixtures_Gernert.pdf)).
This thesis is not yet a curated source record. The relevant pages were checked.

Bose and Cole (`bose_1970`, doi:10.1063/1.1672658), p. 141, report the fitted
$B_{12} = -30.1$ cm³/mol. Their estimate that the derived mixture pressure virials
are uncertain by roughly 2–3 cm³/mol is not a separately stated standard uncertainty
of that cross coefficient. Its fit fixes the pure coefficients, so their uncertainties
and correlations matter. The original page was visually checked. Gernert's Table A2.1,
p. 198, lists five Bose and Cole cross-virial data points with zero used in the fit;
the table image was checked. The comparison is therefore useful as an independent
validator, even though its proposed observable and tolerance need revision.

**Requested resolution:** define comparison at the measured compositions, with the
uncertainty of the pure coefficients and of the cross-coefficient inference included.
Establish the acceptance rule before implementation and include the model's stated
limitations. A disagreement must be reportable; it must not force a tolerance chosen
after observing the residual.
**Status:** open, Q11.

(oep4-review-r8)=
#### R8 (P2) — Complex optics conflicts with the promised pure refractivity limit

**Where:** §4.1's inversion, §6's complex susceptibility and pure limit, and §8's
roundoff-level pure-limit check.

The proposed $\chi$ is complex and §4.1 obtains $n-1$ by inverting it. Existing
`Model.refractivity` instead evaluates the real-valued law from the real
polarizability. In an absorbing gas that differs even from the real part of the full
complex index. `Model.refractive_index` evaluates that complex index.

For the recommended O3 model at 255 nm, 1 bar and 294 K, the review obtained:

| Quantity | Value |
|---|---:|
| Existing `Model.refractivity` | $4.5393908\times10^{-4}$ |
| Real part of `Model.refractive_index` minus one | $4.5388635\times10^{-4}$ |
| Difference, complex-index real part minus existing refractivity | $-5.2724\times10^{-8}$ |
| Imaginary part of the complex index | $5.6189256\times10^{-4}$ |

O3 uses the existing ideal-density convention because it has no pure-fluid EOS. The
example is a contract check in an absorbing band, not a recommended operating state.
The discrepancy is well above roundoff. OEP 1's complex polarizability derivative
does not remove the distinction between the two refractivity conventions.

**Requested resolution:** define whether mixture `refractivity` preserves the
real-polarizability convention or returns the complex index increment. Specify the
matching pure reference, return type and absorption warnings. Add an explicitly
absorbing constituent to the value and derivative pure-limit checks.
**Status:** open, Q12.

(oep4-review-comments)=
#### Further comments C1–C4

- **C1 — Density derivative units (§4.1, §6).** The displayed derivative is with
  respect to molar density $\rho_s$, while the public method takes number density
  $N_s$. State explicitly that
  $\partial\chi/\partial N_s = (1/N_A)\,\partial\chi/\partial\rho_s$, with units m³.
  The displayed molar derivative is algebraically correct; this comment asks for the
  conversion in the interface contract. Open, Q13.
- **C2 — Arrays, species and zero densities (§6, §8).** Specify how wavelength and
  density arrays broadcast, whether scalar densities can mix with arrays, and how
  missing, extra, aliased or repeated species keys are treated. State what happens
  when a constituent is present in the model but has zero density, including range
  warnings, and when every density is zero so that composition cannot be recovered
  by division. Keep the zero-density optical limit distinct from OEP 1's refusal of
  vacuum in the EOS inverse. Open, Q14.
- **C3 — Static versus optical cross coefficients (§5, §7.3).** Require temperature,
  wavelength/static designation and the physical contribution needed to distinguish
  a molecular static dielectric coefficient from an optical refractivity coefficient.
  The source discussions in §3 correctly make this distinction, but the new quantity
  and selection rules must preserve it. A static molecular value must not become an
  optical cross term merely because it is the only curated value. Open, Q15.
- **C4 — Independently mergeable documentation (§12).** Phases 1 and 2 deliver public
  functionality and therefore need their relevant manual documentation when they
  merge, under OEP 0 and the repository's definition of done. Phase 3 can complete the
  example and combined sensitivity discussion; it cannot postpone all manual coverage
  of the earlier phases. Open, Q16.

#### Checks that support the proposal

The ordered-pair counting, single Lorentz–Lorenz inversion, algebraic derivatives and
leading partial-pressure total-density expansion are consistent. This is a check of
the stated algebra and low-density limit, not qualification of every approximation
over the intended operating domain. Garberoglio and Harvey [5] were consulted for the
virial definitions, frequency dependence and classical-limit formulation.

The review reproduced §3's relative total-density differences at 294 K, using
$N = p/(Z k_B T)$ throughout, consistent with optmat's convention rather than using
the backend's molar density with its historical gas constant:

| Mixture | 10 bar: $N_{\rm pp}/N_{\rm mf}-1$ | 30 bar: $N_{\rm pp}/N_{\rm mf}-1$ |
|---|---:|---:|
| He/Ar 50/50 | $+1.6030\times10^{-3}$ | $+5.5212\times10^{-3}$ |
| Ne/Ar 50/50 | $+3.3684\times10^{-4}$ | $+1.6701\times10^{-3}$ |
| N2/O2/Ar air | $-1.5619\times10^{-3}$ | $-4.0824\times10^{-3}$ |

For this table, air uses the CIPM major fractions `(0.780848, 0.209390, 0.009332)`
renormalised to sum to one. `mf` means the `HEOS` multi-fluid mixture;
$N_{\rm pp}$ is the sum of pure-fluid densities at $x_s p$.

The separate air validation probe compared the multi-fluid $Z$ with the dry-air
expression in `picard_2008`'s `appendix_a` formula on a 13-by-6 grid: temperatures
288.15–300.15 K in 1 K steps and pressures 60–110 kPa in 10 kPa steps. It used
solprop's then-current three-component fractions `(0.7808, 0.2095, 0.0093) / 0.9996`,
with zero water vapour in the CIPM expression. The largest absolute difference was
$8.57065\times10^{-6}$, at 288.15 K and 110 kPa, below the proposed
$1.5\times10^{-5}$ criterion. This is an encouraging sampled comparison, not a proof
over the continuous domain or a completed uncertainty budget for composition changes.
The pressure and temperature endpoints were included. R6 records the 39-pair support
count reproduced by the backend parameter queries.

#### Primary-source verification and limitations

- **Freshly checked:** Achtermann et al.'s Table II image and its units on p. 5675;
  Bose and Cole's p. 141 cross-coefficient fit and uncertainty discussion; Gernert's
  virial mixing discussion and Table A2.1 image; the relevant passages of Tkaczuk
  et al. [4] and Garberoglio and Harvey [5]. The two NIST-hosted articles are publicly
  available at [Tkaczuk et al.](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929167)
  and [Garberoglio and Harvey](https://nvlpubs.nist.gov/nistpubs/jres/125/jres.125.022.pdf).
- **Checked against repository records, not freshly against all original pages:**
  the selected like coefficients and their ratios, the Hourri CO2–SF6 comparison,
  and the CIPM formula and composition. The backend capability and phase-contract
  checks also used the official
  [CoolProp mixture documentation](https://coolprop.org/fluid_properties/Mixtures.html).
- **Not verified in full:** Buckingham and Pople [1], the Landolt–Börnstein volume
  [2, 7], and every underlying EOS source and fitted dataset. No new optical
  unlike-pair measurements were acquired or curated. No unlike-term uncertainty
  distribution, temperature extrapolation, GVD sensitivity budget or full continuous
  domain was established by this review. The Chebyshev degree claims of §4.1 were
  not independently reproduced.
- **Initial review checks:** `tests/test_proposals.py` and `tests/test_register.py`
  passed, 47 passed and one skipped; `git diff --check` was clean. The initial
  spelling-check attempt could not run because `codespell` was absent from the
  environment. No full gate, extraction run, refit or documentation build was claimed
  for that read-only review. Validation of this written record is reported separately
  below.

**Handoff:** resolve or explicitly defer Q5–Q16 under OEP 0, recording the resolution
beside each question and preserving this dated review as evidence of what was checked.
The review does not approve a replacement combining rule, an uncertainty model, a
truncation cutoff or an EOS fallback. Those choices still need the maintainer's
decision, supporting evidence and, when implemented, the required register entries.

#### Validation of the written review record (2026-09-30)

The proposal-only checks passed after recording the review:

- `git diff --check`;
- `uv run codespell docs/oeps/oep-0004-mixture-index.md`;
- `uv run pytest tests/test_proposals.py tests/test_register.py`: 47 passed, one
  skipped (the process OEP has a subject-specific structure);
- `uv run poe docs`: the Sphinx build with `-W -n --keep-going`, including the
  gallery examples, succeeded.

The existing locked development dependencies were installed to run these checks,
resolving the initial missing `codespell`. No package dependencies, data, model
coefficients or fitted results were changed. The full gate, library extraction and
fit checks are not required for this proposal-only record and were not run.
- 2026-10-01: The author settled Q17-Q19 (TC-792 to TC-796). Every phase merged and
  reported; status Final, approved by the author, and the OEP archived.
