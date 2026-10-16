# Dispersion models

A model gives the per-particle polarizability $\alpha$ of one gas as a function of
vacuum wavelength (`optmat.gas.models.base.Model.polarizability`); its `refractivity` gives
$n - 1$ at any state through Lorentz-Lorenz. Every model states the wavelength range
its source supports and warns (`ModelRangeWarning`) outside it.

For a first calculation, start with the {doc}`quickstart`. This chapter follows
the shared {doc}`conventions`; {doc}`user-guide/linear` gives the reading order
through model construction, accuracy, absorption, and mixtures.

## Published formulas

A published formula is recorded in its source's record as a `[[formulas]]` entry: the
equation as printed, its variable and units, the state it refers to, and a `check`
saying whether it reproduces the values the paper itself prints. Entries with a `form`
are models (`optmat.gas.models.formulas.published_models`). Their coefficients are kept **as
printed**; a misprinted coefficient is corrected by a `corrections` entry with the
printed and corrected values, the reason and the evidence, like a dataset erratum:

```toml
[[formulas.corrections]]
term = "C"
index = 1
printed = "50.578"
corrected = "59.578"
reason = "..."
evidence = "..."
```

| Form | Formula ($\lambda$ vacuum wavelength, $s = 1/\lambda$ in µm$^{-1}$, $\nu$ in cm$^{-1}$) | Terms |
|---|---|---|
| `dn_sigma2` | $n - 1 = \mathrm{scale}\,(A + \sum_i B_i/(C_i - s^2))$ | `scale`, `A`, `B`, `C` |
| `dn_wavenumber2` | $n - 1 = \mathrm{scale}\,(A + \sum_i B_i/(C_i - \nu^2))$ | `scale`, `A`, `B`, `C` |
| `n2_sellmeier` | $n^2 - 1 = A + \sum_i B_i\lambda^2/(\lambda^2 - C_i)$, $\lambda$ in µm | `A`, `B`, `C` |
| `alpha_series` | $\alpha = \sum_k c_k\,\omega^{2k}$ in atomic units | `c` |
| `molar_refractivity_series` | $A_R = \sum_k a_k\,\lambda^{-2k}$ in cm$^3$/mol, $\lambda$ in nm | `a` |
| `oscillator` | $\alpha = \sum_i f_i/(E_i^2 - \omega^2)$ in atomic units; the entry's `energy_unit` is `hartree` (default), `eV` or `rydberg` | `f`, `E` |

A refractivity formula's `state` and `density` say what the source meant: the
polarizability is recovered with inverse Lorentz-Lorenz at that state's real-gas (or
ideal-gas) density, so formulas referred to different states compare directly.

`kind` distinguishes formulas from papers (`published`) from unpublished formulas
kept for comparison (`legacy`, source `luna_physdata`). The latter reproduce a
recorded implementation snapshot; they should not be read as statements about the
current behaviour of a propagation package. The source record identifies their
origin and coefficients.

## Frequency derivative

Every model gives the derivative of its complex polarizability with frequency,
`dpolarizability_domega(wavelength)`, in C m² V⁻¹ s at fixed state: the quantity a
propagation code needs for the group index and its frame derivatives (OEP 1, section
4.1). Each class gives it analytically:

- an oscillator model, the derivative of its sum of poles, with helium's retardation
  term and, for ozone, the measured absorption's derivative branch by branch: the
  spline's own inside the grid, the differentiated Kramers-Kronig sum above it, and
  the interpolated cross section's slope (the mean of the two at a knot);
- a published formula or a refractiveindex.info entry, the formula's derivative in a
  variable proportional to frequency (1/λ), times the Jacobian of inverse
  Lorentz-Lorenz at its reference density.

It is 0 in the static limit (NaN for the refractiveindex.info entries that are NaN
there) and warns outside the model's range as `complex_polarizability` does.

```python
import math
from optmat.gas.models.registry import model

argon = model("fit/Ar")
assert argon.dpolarizability_domega(math.inf) == 0
slope = argon.dpolarizability_domega(800e-9)
assert slope.real > 0 and slope.imag == 0
```

## Comparison: refractiveindex.info

The gas entries of the refractiveindex.info database are kept, from a pinned snapshot,
as a comparison source only (`optmat.gas.comparison.database_models`,
`optmat.gas.comparison.tabulated_entries`). Its Bideau-Mehu formulas agree with the curated ones
to a few 10$^{-7}$ (it stores the products $1.2055\times10^{-2} f_i$ rounded), and it
references Ermolov's helium formula to 1 atm where the paper says 1 bar.

## Choosing a model

`optmat.gas.models.registry.models(gas)` lists every model of a gas: the fit (`fit/<gas>`), for
the infrared-active molecules its compact form (`compact/<gas>`, a few effective
infrared poles in place of the bands; {doc}`fitting`), the published formulas
(`<source>/<formula>`) and the refractiveindex.info entries.
`optmat.gas.models.registry.default_model(gas)` returns the recommended one, chosen with
its reason in `data/gas/defaults.toml`: today the fitted model for every gas
({doc}`gas_models` says what each rests on). `optmat.gas.models.registry.refractivity(gas,
wavelength, state)` is the one-call entry point, and
`optmat.gas.models.registry.refractive_index` gives the complex index (for ozone, whose
model is absorbing):

```python
from optmat.conditions import State
from optmat.gas.models.registry import refractivity

refractivity("He", 800e-9, State(30e5, 293.15))  # n - 1 at 30 bar, 20 °C
```

## Mixtures: dry air

A mixture is not a material: its make-up is curated as the mole fraction of each
constituent (TC-543), and {func}`optmat.gas.composition.composition` returns it, the
most abundant constituent first, each fraction with its source dataset. The
recommended dry air is the CIPM-2007 composition (Picard et al., Metrologia 45, 149
(2008); TC-548):

| Constituent | Mole fraction | Source |
|---|---|---|
| N2 | 0.780848 | CIPM-2007 |
| O2 | 0.209390 | CIPM-2007 |
| Ar | 0.009332 (u = 3e-6) | CIPM-2007 |
| CO2 | 0.00040 | CIPM-2007 (a conventional value) |
| Ne, He, CH4, Kr, H2, N2O, CO, Xe | 18.2, 5.2, 1.5, 1.1, 0.5, 0.3, 0.2, 0.1 × 10⁻⁶ | US Standard Atmosphere 1976, as CIPM-2007 adopts it |

The fractions sum to 1 − 2.9e-6, the table's rounding. A calculation using only
some constituents must specify which are retained and how their fractions are
renormalised.

```python
from optmat.gas.composition import composition

for constituent in composition("dry air"):
    print(constituent.material, constituent.mole_fraction, constituent.source)
```

Two reduced compositions retain and renormalise CIPM-2007 constituents:
`"dry air (N2, O2, Ar)"` (99.957 % of the molecules; the air of the Raman response)
and `"dry air (N2, O2)"`. `composition("dry air")` returns the full table as printed.
The older conventional oxygen fraction 0.2095 differs from CIPM-2007's 0.20939;
choose a composition explicitly when reproducing older calculations (TC-759).

**The air model.** `optmat.gas.air.air_model(name)` is a named composition as a
dispersion model: its polarizability per molecule is the mole-fraction sum of the
constituents' recommended models, $\alpha_\text{air} = \sum_s x_s \alpha_s$, with the
derivative likewise, the intersection of their intervals (widened where a trace
constituent allows, below), and an identity bound to theirs. The default is `"dry air (N2, O2, Ar)"`; `"dry air (CIPM-2007)"` holds all
twelve constituents, renormalised by their printed sum. Its density at a state is
CoolProp's pseudo-pure dry air (Lemmon et al. 2000, `lemmon_2000`), which holds the
unlike-pair interactions: the constituents' partial-pressure densities are
1.9 × 10⁻⁴ below it at 1 bar and 294 K, 1.8 × 10⁻³ at 10 bar and 6.0 × 10⁻³ at 50 bar.
This air model does not include unlike-pair refractivity terms. For a labelled
mixture with pair terms and its own density law, use {doc}`mixtures`.

| Composition | Against the full composition (polarizability) | Against Ciddor (1996), n − 1 at 15 °C, 1 atm, 400 umol/mol CO2 |
|---|---|---|
| `dry air (CIPM-2007)` | — | −2.7 × 10⁻⁸ at 320 nm; within ±1.1 × 10⁻⁸ from 400 to 1400 nm |
| `dry air (N2, O2, Ar)` | −1.88 to −1.95 × 10⁻⁴ (the CO2 and traces left out), 320–1400 nm | −8.4 to −4.1 × 10⁻⁸ |
| `dry air (N2, O2)` | +1.50 to +1.86 × 10⁻⁴ (argon and the rest left out), 320–1400 nm | +2.6 to +5.2 × 10⁻⁸ |

Since OEP 3's refits (TC-797) the intersection of the full composition's constituents'
intervals is 320-1313 nm, with gaps at 1123-1137, 1149-1177 and 1187.4-1187.9 nm, where
its trace CH4 absorbs as a pure gas. A trace constituent can leave the mixture
transparent where its own model is not qualified, so the intersection widens where the
mixture passes OEP 3's checks at its mole fractions (Q13): its absorption below
0.1 m⁻¹, no pole of an unqualified constituent, and their representation errors, from
each molecule's gap table, within 0.1 of the mixture's statistical sigma in the
polarizability and 1 % of its $D_2$ (TC-809). The full composition gains 1499-1516 and
1548-1562 nm; its trace gases' other gaps stay closed, because their band-rule poles
fill them. Below N2O's 320 nm, whose ultraviolet absorption is not assessed, nothing is
widened. The reduced compositions keep one interval.

Ciddor's standard dry air is an independent refractometric source (its stated
uncertainty 2–5 × 10⁻⁸ in n); the full composition agrees with it within that. The
reduced compositions' difference from it is their omission (TC-759).

```python
from optmat.conditions import State
from optmat.gas.air import air_model

air = air_model()  # "dry air (N2, O2, Ar)"
n_minus_1 = air.refractivity(800e-9, State(1e5, 294.0))
assert 2.6e-4 < n_minus_1 < 2.8e-4
```
