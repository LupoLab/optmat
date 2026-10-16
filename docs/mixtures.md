# Gas mixtures at finite density

This chapter extends the pure-gas calculation in the {doc}`quickstart`. Read
{doc}`conventions` and {doc}`models` first, especially the distinction between
polarizability and refractive index. The API is {mod}`optmat.gas.mixture` with
{mod}`optmat.gas.mixture_eos`; the {doc}`review register <register>` records the
pair-model decisions.

A gas-filled fibre often holds a mixture: helium and argon, neon and argon, air, or a
molecular gas in a noble-gas buffer. At the tens of bar of such fills the refractive
index of a mixture is not the sum of its constituents' indices. It is the Lorentz-Lorenz
law of the whole mixture, including the pair terms between unlike molecules, at the
mixture's own density. This page covers both
([OEP 4](oeps/archive/oep-0004-mixture-index.md)): the optics at stated densities,
and the equation of state that gives the densities at a pressure and temperature.

```python
import numpy as np
from optmat.gas.mixture import mixture_optics

he_ar = mixture_optics(["He", "Ar"])
n = 1.25e26  # m^-3 of each: about 10 bar together at 294 K
dn = he_ar.refractivity(800e-9, {"He": n, "Ar": n})  # n - 1
chi = he_ar.susceptibility(800e-9, {"He": n, "Ar": n})  # complex n^2 - 1
dchi = he_ar.dsusceptibility_domega(800e-9, {"He": n, "Ar": n})  # s
per_species = he_ar.dsusceptibility_ddensity(800e-9, {"He": n, "Ar": n})  # m^3
spread = he_ar.refractivity_uncertainty(800e-9, {"He": n, "Ar": n})
assert 1.4e-3 < float(dn) < 1.5e-3 and float(spread) < 2e-4 * float(dn)
assert [p.status for p in he_ar.pairs] == ["curated", "curated", "estimate"]
bound = he_ar.density_bound({"He": 0.5, "Ar": 0.5})  # m^-3, the optical domain
```

With the density from the mixture's equation of state:

```python
import warnings
from optmat.conditions import State
from optmat.gas.mixture import mixture_optics
from optmat.gas.mixture_eos import mixture_equation_of_state

he_ne = mixture_equation_of_state({"He": 0.5, "Ne": 0.5})
assert (he_ne.method, he_ne.status) == ("virial", "validated")
densities = he_ne.partial_densities(State(pressure=1e6, temperature=294.0))  # m^-3
dn = mixture_optics(["He", "Ne"]).refractivity(800e-9, densities)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # a provisional equation warns when it is built
    air = mixture_equation_of_state("dry air (N2, O2, Ar)")
assert (air.method, air.status) == ("multi-fluid", "provisional")
```

## The law

In the low-density cluster expansion the Lorentz-Lorenz function of a mixture is the
sum over its molecules plus the sum over its pairs (Buckingham and Pople 1955):

$$
y = \frac{n^2 - 1}{n^2 + 2} = \sum_s A_s\rho_s + \sum_{s,t} B_{st}\,\rho_s\rho_t,
\qquad A_s = \frac{N_A\alpha_s}{3\varepsilon_0},
$$

with $\rho_s = N_s/N_A$ the molar density of species $s$ and $\alpha_s$ its
polarizability from its recommended model ({doc}`gas_models`). For one species
$B_{ss}$ is that gas's refractivity virial coefficient $B_R$, so the law is optmat's pure
law exactly ({doc}`conventions`). The mixture's $y$ is inverted once,
$\chi = 3y/(1 - y)$ and $n - 1 = \chi/(1 + \sqrt{1 + \chi})$. Summing the constituents'
indices, or their susceptibilities, is wrong at the order of the pair terms.

- **Real and complex.** `refractivity` is in the pure models' real-polarizability
  convention, and for one species equals `Model.refractivity` at that density.
  `susceptibility` uses each constituent's complex polarizability, and for one species
  $\sqrt{1 + \chi}$ equals `Model.refractive_index`. The two differ only for an absorbing
  constituent: for O3 at 255 nm, 1 bar and 294 K, the real part of the complex index lies
  $5.27 \times 10^{-8}$ below the real convention's $n - 1$, which is $(\operatorname{Im}\chi)^2/24$.
  An absorbing constituent warns, as its pure model does.
- **Densities.** `densities` maps each species to its number density in m⁻³. The keys
  must be exactly the mixture's species: a missing or extra species is refused by name.
  The wavelength and the densities broadcast by numpy's rules. Two labelled populations
  of one species are combined by the caller first: they are one species in the pair
  sum.
- **Derivatives** are closed forms. With $\partial\chi/\partial y = 3/(1-y)^2$,
  $\partial_\omega\chi$ takes each constituent's `dpolarizability_domega`, and
  $\partial\chi/\partial N_s = (3/(1-y)^2)\,[A_s + 2\sum_t B_{st}\rho_t]/N_A$.

## The pair terms

Every pair of a mixture has a `PairTerm` (`MixtureOptics.pairs`, like terms first) with
its status, source and reason.

- **Like terms** are each gas's curated $B_R$ (`data/gas/refractivity_virial.toml`),
  dispersing with the gas, $B_{ss} = b_s A_s$. A static (dielectric) value would be taken
  as an optical term only where its induced part is stated; every present choice is
  optical.
- **O2's like term** is Egan et al.'s published estimate, $b = 0.210$ cm³/mol at
  633 nm with 50 % uncertainty: `provisional`, with its uncertainty in the
  sensitivity. Its pairs are estimates, and its mixtures `provisional`.
- **A gas without a curated $B_R$** (D2, N2O, CO, O3) is `unknown`: it enters with
  $b = 0$ and the declared scale $|b| \le 2.5$ cm³/mol, and its mixtures are
  `provisional`.
- **Unlike terms.** No unlike-pair refractivity coefficient is curated: none is known to
  have been measured for the pairs of interest (the search is recorded as TC-775). Each
  is estimated as the mean of its like terms, $B_{st} = \tfrac12(B_{ss} + B_{tt})$, each
  like term dispersing with its own $A$, with status `estimate`, or `unknown` where a
  like term is. Its standard uncertainty is $\max(|B_{ss}|, |B_{tt}|)$, taking an unknown
  like term at its scale. That is a sensitivity assumption, not a statistical
  uncertainty: cancellation between like terms of opposite sign, helium's and argon's,
  cannot shrink it.
- **The sensitivity.** `refractivity_uncertainty` combines the estimated and unknown
  terms' uncertainties linearly, as fully correlated, and gives their effect on $n - 1$.

At 800 nm and 294 K, the ideal-gas density at each pressure, equal mole fractions (air:
CIPM-2007's N2, O2 and Ar, renormalised):

| Mixture | Pressure | $n - 1$ | Its sensitivity | Relative | GVD, fs²/m | Its sensitivity |
|---|---:|---:|---:|---:|---:|---:|
| He/Ar | 10 bar | 1.441e-3 | 2.2e-7 | 1.5e-4 | 96.8 | 0.016 |
| He/Ar | 30 bar | 4.326e-3 | 2.0e-6 | 4.6e-4 | 290.8 | 0.14 |
| Ne/Ar | 10 bar | 1.588e-3 | 2.2e-7 | 1.4e-4 | 101.1 | 0.016 |
| Ne/Ar | 30 bar | 4.767e-3 | 2.0e-6 | 4.2e-4 | 303.6 | 0.14 |
| He/Ne | 10 bar | 4.656e-4 | 1.4e-8 | 3.0e-5 | 14.2 | 0.00042 |
| He/Ne | 30 bar | 1.397e-3 | 1.2e-7 | 8.9e-5 | 42.5 | 0.0037 |
| Air | 10 bar | 2.660e-3 | 8.6e-8 | 3.2e-5 | 205.9 | 0.0063 |
| Air | 30 bar | 7.989e-3 | 7.7e-7 | 9.7e-5 | 619.1 | 0.057 |

Air's sensitivity comes mainly from oxygen's estimated like term and its pairs. The GVD
columns are by differences of `refractivity` and `refractivity_uncertainty` in
frequency.

## The domain

The law is second order in density. `density_bound(composition)` gives the optical
domain, the lower of two densities:

- **the third term.** Where the neglected third term reaches 10⁻³ of $y$:
  $\rho = \sqrt{10^{-3}\,A_R/|C_R|}$, with the largest $|C_R/A_R|$ of the species present.
  $C_R$ is measured for Ar, Kr, Xe, H2, N2, CH4, CO2 and SF6 (Achtermann et al. 1991 and
  1993), with a measured bound for He. Neon, whose $C_R$ is not measured, takes the
  declared scale $|C_R| \le 600\,B_R^2/A_R$. A gas with neither $B_R$ nor $C_R$ gives no
  third-term bound.
- **the estimated and unknown terms.** Where their uncertainty reaches 10⁻³ of $y$.

Both are taken as the smallest over the qualified intervals. The bounds are about 6960
mol/m³ for argon and 6670 for N2 (the third term), 2450 for 50/50 He/Ar
(the pair terms), 6150 for air (oxygen's third term, at the declared scale on its
estimated $B_R$), and 400 for a pure gas with an unknown like term. Beyond the bound the
optics warn (`ModelRangeWarning`). In wavelength, a mixture holds on the intersection
of its constituents' qualified intervals (`MixtureOptics.intervals`), and warns outside
it whatever the densities.

## Identity

A mixture's identity ({doc}`identity`) binds its constituents' models by their digests,
with the pair rules and the third terms. Its status is `provisional` where any term is
`unknown`, and `validated` otherwise. No unlike term has been checked against a
measurement, and the page states that rather than implying it.

## How much each piece matters

For a 50/50 helium-argon fill at 294 K and 800 nm, the relative change in $n - 1$ and in
the GVD when one piece of the law is left out or replaced, against the full law (the
default equation of state and the mixture optics). The gallery example
{ref}`sphx_glr_auto_gallery_accuracy_plot_mixture_he_ar.py` plots it from 1 to 60 bar.

| Pressure | $n - 1$ | Dilute law | No pair terms | Partial-pressure density | Multi-fluid density | Unlike-pair estimate |
|---|---:|---:|---:|---:|---:|---:|
| 10 bar | 1.436e-3 | −1.1e-3 / −2.1e-3 | −1.5e-4 / −1.6e-4 | +6.0e-3 / +6.3e-3 | +2.1e-3 / +2.1e-3 | ±1.5e-4 |
| 30 bar | 4.283e-3 | −3.3e-3 / −6.2e-3 | −4.4e-4 / −4.8e-4 | +1.8e-2 / +1.9e-2 | +5.7e-3 / +5.7e-3 | ±4.5e-4 |
| 60 bar | 8.492e-3 | −6.5e-3 / −1.2e-2 | −8.6e-4 / −9.4e-4 | +3.6e-2 / +3.8e-2 | +9.1e-3 / +9.2e-3 | ±9.0e-4 |

Each cell is the change in $n - 1$ / in the GVD; the last column is the estimate's
declared uncertainty in $n - 1$.

- **The density rule matters most.** The partial-pressure rule leaves out the unlike
  pairs' second virial and, beyond first order, also moves the composition toward argon:
  together 6 × 10⁻³ at 10 bar.
- **The density method is the next.** The two methods differ by 2 × 10⁻³ at 10 bar,
  about ten times the optics' own pair terms and their uncertainty. Since the
  multi-fluid model misses the measured cross virials (TC-779), the default takes the
  fallback, which misses Blancett et al.'s data by up to 7 × 10⁻⁴ in $Z$ at 28 bar. So
  the density, not the optics, sets the accuracy of a helium-argon fill.
- **Lorentz-Lorenz against the dilute law** is 1.1 × 10⁻³ in $n - 1$ at 10 bar and twice
  that in the GVD.
- **The pair terms** are 1.5 × 10⁻⁴ at 10 bar. The declared uncertainty of the unlike
  estimate is about the same size, and it grows linearly with pressure.

## The density

`mixture_equation_of_state(composition, method=None)` gives a mixture's equation of
state at a fixed composition (OEP 4, sections 4.2 and 7.4). `number_density(state)` is
the total density, and `partial_densities(state)` gives each species' $x_s N$, the keys
the optics take. `pressure(number_density, temperature)` is the inverse, and
`compressibility` and `second_virial` complete it. The densities follow the CODATA
convention of {doc}`transport`'s equations of state: $N = p/(Z k_B T)$ with the
equation's own $Z$. `molar_mass` is the composition's, $\sum_s x_s M_s$, kg/mol, with
each constituent's equation-of-state molar mass, whichever the method.

There are three methods, chosen per mixture by the evidence for its pairs:

- **`multi-fluid`**: CoolProp's multi-fluid Helmholtz model, where every pair has a
  binary model fitted to data of that pair. It reduces to each pure fluid's equation at
  $x_s = 1$.
- **`virial`**: each species' own equation at its partial density, plus the measured
  cross second virial coefficients,
  $p = \sum_s p_s(x_s\rho, T) + RT\rho^2\sum_{s\ne t}x_s x_t B_{st}(T)$. It preserves
  the composition exactly and holds while $|\rho\sum_{s\ne t}x_s x_t B_{st}| \le 10^{-2}$,
  beyond which it warns. Each pair's $B_{st}$ is a measured $B_{12}$, or a measured excess
  $E$ with $B_{12} = E + (B_{11} + B_{22})/2$ from the reference equations' pure $B$.
  The method holds only at the temperatures where they are measured, and refuses others.
- **`partial-density`**: the same with every unknown cross coefficient zero, the
  composition-preserving form of the partial-pressure rule. It is `provisional`.

**The comparison.** Each multi-fluid model is compared with every curated measured
cross, excess and mixture second virial, and every mixture density or compressibility
factor, between 250 and 350 K and up to the mixture's optical domain. A point passes
within twice the combined uncertainty; a pair qualifies when at least 95 % of its points
pass in each kind, with at least three points. The rule was fixed before the comparison
(TC-778). The fallback is compared with the densities in the same way, within its
truncation bound (TC-779).

| Pair | Multi-fluid model | Its comparison | Fallback's $B_{12}$ | Fallback against densities |
|---|---|---|---|---|
| He-Ar | Tkaczuk et al. 2020 | $E$ about 10 cm³/mol low; Z up to 1e-2 low at 100 bar: fails | Brewer and Vaughn's $E$ | Blancett et al.'s Z: 16 of 78 pass, up to 7e-4 low: fails |
| Ne-Ar | Tkaczuk et al. 2020 | $E$ 3.4 against 13.0 cm³/mol: fails | Brewer and Vaughn's $E$ | Kestin et al.'s densities: 7 of 8 pass: fails |
| He-Ne | Tkaczuk et al. 2020 | $E$ −2.6 against 0.76 cm³/mol; Z 2e-3 low at 36 bar: fails | Brewer and Vaughn's $E$ | Vogl and Hall's Z: 35 of 36 pass: qualified |
| N2-O2 | Gernert 2013 | $E$ −0.34 against −1.0(2) cm³/mol: fails | Martin et al.'s $E$ | no density data |
| Ar-N2 | Gernert 2013 | $E$ −0.16 against −0.57(10): fails | Brewer and Vaughn's $E$ | no density data |
| Ar-O2 | Gernert 2013 | 3 of 4 pass: fails | Martin et al.'s $E$ | no density data |
| Ar-CH4, CH4-CO2, CH4-N2 | Kunz and Wagner 2012 | Martin et al.'s $E$, 3 of 3 each: qualified | — | — |

Tkaczuk et al. fitted their models to densities and phase equilibria and did not fit the
second virials. CoolProp reproduces their check values (Table 8), so the gap in the
cross virials is their models', not CoolProp's.

**The defaults.**
- **He/Ne:** `virial`, validated.
- **He/Ar and Ne/Ar:** `virial`, provisional, with the deviations above. The fallback
  is five to ten times closer to Blancett et al.'s helium-argon data than the
  multi-fluid model.
- **Dry air:** `multi-fluid`, provisional, at every temperature (the author's decision,
  TC-796). Its $Z$ agrees with CIPM-2007's within 8.6 × 10⁻⁶ over 600–1100 hPa and
  17–27 °C, where CIPM-2007's $u(Z)$ is 1.5 × 10⁻⁵. The fallback agrees within
  1.8 × 10⁻⁶, but holds only at Martin et al.'s 290–320 K; it can be asked for by name
  (`method="virial"`). By the rule of TC-778 alone, the fallback would have been the
  default.
- **The Kunz and Wagner pairs above:** `multi-fluid`, validated.
- **Other mixtures:** the partial-density rule, provisional.

A provisional equation warns when it is built. Any method can be asked for by name where
the pairs support it, and the identity records the method and status.

**Air two ways.** The labelled N2/O2/Ar fill uses this equation of state, and
`Gas.Air`'s pseudo-pure dry air (Lemmon et al. 2000) serves the air models of
{doc}`gas_models`. At 10 bar and 294 K the multi-fluid N2/O2/Ar density lies 2.0 × 10⁻⁴
below the pseudo-pure one, and the fallback 1.8 × 10⁻⁴ above the multi-fluid one. Their
molar masses differ too: 28.9597 g/mol for the N2/O2/Ar fill, against CoolProp's
28.96546 g/mol for the pseudo-pure air (TC-832).

**Phases.** The multi-fluid method checks the mixture's own phase with CoolProp's flash,
and a state that is not a single gas phase raises `StateError`. The fallback checks each
constituent at its partial density against its own saturated vapour. That is a necessary
condition, not a sufficient one: its mixtures are not checked for phase separation.

This page documents what [OEP 4](oeps/archive/oep-0004-mixture-index.md) proposed and
implemented: its theory, decisions and phase reports are the record.
