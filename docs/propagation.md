# Using optmat from a propagation code

A pulse-propagation code can obtain gases' and solids' linear
and nonlinear response, densities and transport from optmat. It can also record
which models it used. This page gathers what such a code calls, in the units it gets them,
and what it should store
({doc}`OEP 1 <oeps/archive/oep-0001-propagation-interface>`). Each snippet runs in the test suite.

## The surface

The supported names are each module's `__all__`, and a test keeps a snapshot of their
signatures (`tests/data/public_api.json`). What a propagation code needs:

| Need | Call |
|---|---|
| A gas's linear response | `optmat.gas.models.registry.default_model(gas)`: `polarizability`, `complex_polarizability`, `dpolarizability_domega`, `refractivity(wavelength, state)`, `refractive_index(wavelength, state)` |
| Air | `optmat.gas.air.air_model()` |
| A gas's density | `optmat.gas.eos.equation_of_state(gas)`: `number_density(state)`, `pressure(N, T)`, `phase_domain(T)`, `speed_of_sound(state)`, `molar_mass`; air's: `optmat.gas.eos.air_equation_of_state()` |
| A mixture | `optmat.gas.mixture.mixture_optics(species)` and `optmat.gas.mixture_eos.mixture_equation_of_state(composition)` ({doc}`mixtures`) |
| The Kerr effect of a gas | `optmat.gas.nonlinear.hyperpolarizability_models()[gas]`: `gamma3`, `n2`, `gamma_si`; where its fit is supported and the gas's resonance limits: `support` |
| The Raman response of a gas | `optmat.gas.raman_response.raman_response(gas, temperature)`: `response`, `spectrum`, `groups()`; its status alone, or None for none: `raman_response_status(gas)` |
| Ionisation | `optmat.gas.ionisation_rates.ionisation_rate_model(gas, wavelength)`, and the molecules' `molecular_rate_model` |
| Two-photon ionisation of a noble gas | `optmat.gas.two_photon.two_photon_ionisation_model(gas, polarization)`: `cross_section`, `available`, `region`, `rate`, `loss_coefficient`, `pair_cross_section`; NaN where unavailable, never an error on a grid ({doc}`two_photon`) |
| Transport | `optmat.gas.transport.transport_model(gas)`: `viscosity`, `heat_capacity`, `molar_mass` |
| Fused silica | `optmat.solid.bundle.silica_bundle()`, and the solids' `optmat.solid.models.registry` |

## Units

Everything is SI, with these conventions:
- **Wavelengths** are vacuum wavelengths in metres; `inf` is the static limit.
- **Linear response** is per particle: a polarizability in C m²/V. Its frequency
  derivative is `dpolarizability_domega`, in C m² s/V.
- **Densities** are number densities in m⁻³, and pressures are in Pa.
- **States** are a `State(pressure, temperature)`. A refractivity at a state takes the
  real-gas density by default (`ideal=True` for the ideal gas).

```python
import numpy as np
from optmat.conditions import State
from optmat.gas.eos import equation_of_state
from optmat.gas.models.registry import default_model

argon = default_model("Ar")
wavelengths = np.array([400e-9, 800e-9, 1.6e-6])
alpha = argon.polarizability(wavelengths)  # C m^2/V, per atom
dalpha = argon.dpolarizability_domega(wavelengths)  # C m^2 s/V
state = State(pressure=1e6, temperature=294.0)
density = equation_of_state("Ar").number_density(state)  # m^-3
dn = argon.refractivity(wavelengths, state)  # n - 1, Lorentz-Lorenz with B_R
assert alpha.shape == dalpha.shape == dn.shape == (3,)
assert 2.4e26 < density < 2.6e26
```

**The third order.**
- `gamma3` is the per-molecule coefficient $\gamma_3 = \gamma_T/(6\varepsilon_0)$
  in m⁵/V², with $P = \varepsilon_0\gamma_3 N E^3$ for the real field.
- `n2` is the dilute-limit self Kerr index of linearly polarised light at a state, in
  m²/W.
- `gamma_si` is the molecular hyperpolarizability itself.

The factorial in this definition is essential when converting between `gamma3`,
the Taylor-series hyperpolarizability and the susceptibility. See
{doc}`user-guide/nonlinear-conventions`.

```python
from optmat.gas.nonlinear import hyperpolarizability_models

kerr = hyperpolarizability_models()["Ar"]
gamma3 = kerr.gamma3("kerr", 800e-9)  # m^5/V^2
n2 = kerr.n2(800e-9, State(pressure=1e5, temperature=294.0))  # m^2/W
assert 1e-24 < float(n2) < 1e-22
```

## What to store

Store the `identity` of every model a run used ({doc}`identity`). It is the model's
name, the quantity, the material, its sources, its status and a SHA-256 **digest** of
everything its numbers depend on: its parameters, its intervals, the models it is bound
to, and the evaluation code's versions. A later run reproduces the same numbers when the
digests match. A model built from a caller's parameters has the status `user`, and
records what changed and the curated model it started from.

```python
identity = argon.identity
record = {
    "name": identity.name,
    "digest": identity.digest,
    "status": identity.status,
    "optmat": identity.optmat_version,
}
assert identity.kind == "gas_refractive" and len(record["digest"]) == 64
```

**The status.** `validated` is a default with a stated uncertainty. `provisional` and
`untested` say why in `status_reason`. `comparison` is a curated alternative. A code
should report the status of every model it used, and treat `provisional` and `untested`
results as such.

## Where a model holds

- **Intervals.** A dispersion model holds on its qualified `intervals` of vacuum
  wavelength. Outside every one it warns (`ModelRangeWarning`), and between two of
  them the warning names the gap's reason.
- **Unphysical wavelengths.** Every evaluation refuses a wavelength at or below zero
  (`ArgumentError`), and a solid also refuses `inf`; a gas and air take `inf` as the
  static limit. A code whose frequency grid reaches $\omega \le 0$ masks those points
  before it asks for an index. NaN passes through as a missing value (TC-1091).
- **Absorption.** An absorbing gas warns (`ModelAbsorptionWarning`) where its real part
  leaves the absorption out; the complex methods keep it.
- **Phase.** The equations of state refuse a condensed state (`StateError`).
- **The third order.** A hyperpolarizability model warns beyond 1.1 times its largest
  fitted $W_2$. A code that evaluates every four-wave quartet leg by leg reads
  `support`, a record of facts with its own digest: the largest fitted $W_2$, single
  leg and pair sum, with the points that set them, and the gas's resonance limits $E_1$,
  $E_2$ and the ionisation energy, each with its datasets, kind and uncertainty. Two
  derived quantities come with them, `leg_max_in_w2` and `band_edge()`. What to do
  beyond the data is the code's own decision; store the record's digest beside the
  model's ({doc}`user-guide/nonlinear-support`;
  {doc}`OEP 5 <oeps/archive/oep-0005-hyperpolarizability-support>`).
- **Warnings are part of the answer.** A propagation code should surface them, not
  silence them.

```python
import warnings
from optmat.exceptions import ArgumentError, ModelRangeWarning, StateError

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    argon.polarizability(50e-9)  # below argon's qualified range
assert any(issubclass(w.category, ModelRangeWarning) for w in caught)
try:
    equation_of_state("SF6").number_density(State(pressure=3e6, temperature=294.0))
except StateError:
    pass  # SF6 is a liquid at 30 bar and 294 K
else:
    raise AssertionError("a condensed state must be refused")
try:
    argon.polarizability(-800e-9)
except ArgumentError:
    pass  # a negative wavelength is not a frequency a model holds at
else:
    raise AssertionError("an unphysical wavelength must be refused")

support = kerr.support  # facts in m^-1 (W2 in m^-2), with their sources
assert support.model_digest == kerr.identity.digest != support.digest
assert support.leg_max.wavenumber <= support.leg_max_in_w2  # the largest fitted leg
# For any polarisation the bound is lower in wavenumber: 192.1 nm against 189.4 nm.
assert support.band_edge("any") < support.band_edge() < support.e1.wavenumber
```

## The Raman response

`raman_response(gas, temperature)` gives a gas's (or a mixture's) Raman response.
- `groups()` collects its lines into groups that share one decay rate.
- `decay_rate(N)` is defined down to zero density, the collisionless limit.
- Each line has a `line_id` (`"N2/rotation/S/0"`).
- `check_densities(N)` refuses densities outside the width law's range.

A code that integrates the response as damped oscillators can carry one decay per group
({doc}`raman_response`).

`raman_response_status(gas)` reads the status the response would have, without
computing it. It returns None for a gas with no Raman-active constituent.

```python
from optmat.gas.raman_response import raman_response, raman_response_status

assert raman_response_status("Ar") is None
assert raman_response_status("dry air") == ("validated", "")
nitrogen = raman_response("N2", 295.0)
groups = nitrogen.groups()
assert nitrogen.line_id[0] == "N2/rotation/S/0"
assert all(np.all(np.isfinite(g.decay_rate(0.0))) for g in groups)
h = nitrogen.response(np.array([0.0, 1e-13, 1e-12]), 2.5e25)
assert h.shape == (3,)
```

## Ionisation

`ionisation_rate_model(gas, wavelength)` gives the default rate:
- `instantaneous_rate(E)` in 1/s, for the field's instantaneous magnitude in V/m;
- `cycle_averaged_rate`;
- `channel_rates(E)`, the rate of each ionic channel;
- `energy_per_event(E)`, in J: the energy the field gives up per ionisation event, of
  the kind `energy_kind` names;
- `energy_loss_rate(E)` in W per particle.

These rates are for linear polarisation ({doc}`ionisation_rates`).

```python
from optmat.gas.ionisation_rates import ionisation_rate_model

rate = ionisation_rate_model("Ar", 800e-9)
field = 2e10  # V/m
w = rate.instantaneous_rate(field)  # 1/s
energy = rate.energy_per_event(field)  # J
assert rate.energy_kind == "adiabatic"
assert np.isclose(rate.energy_loss_rate(field), w * energy, rtol=1e-12)
```

## Fused silica

`silica_bundle()` pairs the Kerr model and the Raman response that belong together. It
gives the total, electronic and Raman $n_2$, the Raman fraction, and the response $h_R$.
Its `gain_check` compares its Raman gain with a measured one ({doc}`silica`).

```python
from optmat.solid.bundle import silica_bundle

silica = silica_bundle()
n2_total = silica.n2(1030e-9)  # m^2/W
assert np.isclose(
    n2_total, silica.n2_electronic(1030e-9) + silica.n2_raman(1030e-9), rtol=1e-12
)
assert 0.1 < silica.raman_fraction < 0.3
```

## Overrides

`with_parameters(...)` and the classes' `from_parameters` build a model from a caller's
own numbers, checked as a curated one is. The result has the status `user`
({doc}`identity`).

```python
narrow = nitrogen.with_parameters(width_scale=0.8)
assert narrow.identity.status == "user"
assert narrow.identity.base == nitrogen.identity.digest
```

## Transport and mixtures

`transport_model(gas)` gives the dilute viscosity, the ideal-gas $c_p$ and the molar
mass ({doc}`transport`). A mixture's optics and density take its species' number
densities, which `partial_densities` returns ({doc}`mixtures`).

```python
from optmat.gas.mixture import mixture_optics
from optmat.gas.mixture_eos import mixture_equation_of_state
from optmat.gas.transport import transport_model

eta = transport_model("He").viscosity(294.0)  # Pa s
he_ne = mixture_equation_of_state({"He": 0.5, "Ne": 0.5})
densities = he_ne.partial_densities(State(pressure=1e6, temperature=294.0))
chi = mixture_optics(["He", "Ne"]).susceptibility(800e-9, densities)
assert he_ne.status == "validated" and float(chi.real) > 0
```

This page documents what [OEP 1](oeps/archive/oep-0001-propagation-interface.md)
proposed and implemented: its theory, decisions and phase reports are the record.
