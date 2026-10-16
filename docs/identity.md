# Identity, status and validity

A code that stores results computed with optmat needs to know, later, which model
produced them, whether that model has changed since, and how far it was to be trusted.
Every model, response and bundle of the public surface answers all three through its
`identity`, its status, and its qualified intervals ({doc}`OEP 1 <oeps/archive/oep-0001-propagation-interface>`,
sections 4.1, 4.8 and 6.2).

The {doc}`quickstart` shows how to select a model. Read this chapter to decide what
to store with a calculation, then continue through {doc}`user-guide/linear` or
use {doc}`propagation` for the integration interfaces.

```python
from optmat.gas.models.registry import model

n2 = model("fit/N2")
identity = n2.identity
assert identity.name == "fit/N2"
assert identity.kind == "gas_refractive"
assert identity.status == "validated"
assert len(identity.digest) == 64
```

## The identity

`identity` is a {class}`~optmat.identity.ModelIdentity`: the model's name, the quantity
it gives (`kind`), the material, the source records it rests on, its status, and a
**digest**. The digest is the SHA-256 of a canonical record of the model:

- the class, and the `EVALUATION_VERSION` of the class and of every class and shared
  evaluator module its evaluation goes through; a version is raised whenever that code
  changes its output;
- every field of the model: the parameters that fix a state (temperature, wavelength,
  composition) and the values evaluation reads (a fit's terms and offsets, a
  spectrum's grid), not the files that produced them;
- the qualified intervals (below);
- the values a model resolves once from another component, and the digests of the
  components it is bound to: the refractivity virial ratio and a formula's reference
  density, the Kerr fit and linear model behind silica's Raman kernel, the nitrogen
  model behind a solid formula relative to nitrogen.

The record is canonical JSON: sorted keys, no insignificant whitespace, floats in
Python's shortest round-trip form, non-finite floats as `"inf"`, `"-inf"` and `"nan"`,
arrays as lists and a complex number as `[re, im]`. The same model gives the same digest
in every process and on every IEEE-754 platform whose arithmetic produced the same
values.

The digest changes whenever the model's numbers can change. The converse is not
promised: an edit that leaves every number unchanged may still change it. The package
version is recorded beside the digest (`optmat_version`), not in it, so a digest
survives a release that did not touch the model. A code that cannot accept the risk of
an evaluation change escaping the version bump also compares the version.

The test suite checks the versions at sampled arguments using recorded model identities
and golden values. For gas Raman responses, derived widths can differ in their last
bits between numerical libraries. The tests compare their canonical identity records
with an arithmetic tolerance, while checking versions and discrete metadata exactly;
the runtime digest remains an exact hash. Other models retain exact digest checks.
All sampled outputs are compared independently. See {doc}`developers/testing` for
the tolerances and fixture workflow. These checks cannot establish equivalence at
every argument.

`identity.kind` names the quantity (`gas_refractive`, `gas_hyperpolarizability`,
`ionisation_rate`, `raman_response`, `solid`, `kerr`). It is not a dispersion model's
own `kind`, which names its provenance (`published`, `fitted`, ...).

**Records about a model.** A hyperpolarizability model's `support`, the record of where
its fit is supported and of the gas's resonance limits ({doc}`nonlinear`), is not a field
of the model and does not enter its digest: it changes no evaluated number. It is
identified in its own right (`identity.kind` is `gas_hyperpolarizability_support`), and
its digest covers every fact, witness, kind, uncertainty, missing reason and the
digests of what it was derived from: the selection, the fit's configuration and each
selected dataset. Its `model_digest` names the model it describes. A code stores both
digests.

## The status

| Status | Meaning |
|---|---|
| `validated` | A default, with a stated uncertainty |
| `provisional` | A known stopgap, with its reason: CO2's and N2O's Raman responses (TC-542) |
| `untested` | No reference to test it against: SF6's molecular ionisation rate (TC-678) |
| `comparison` | A curated model kept to compare with the default: the other published formulas, the refractiveindex.info entries, the literature's rate forms, the other silica shapes |
| `user` | Built from parameters a caller supplied; a curated model's stated uncertainty does not carry over |

`status_reason` says why, where that needs saying. A mixture's Raman response takes the
least trusted status of its Raman-active constituents.

```python
from optmat.gas.molecular_ionisation_rates import molecular_rate_model
from optmat.gas.raman_response import raman_response

assert molecular_rate_model("SF6", 800e-9).status == "untested"
assert raman_response("CO2", 295.0).status == "provisional"
assert raman_response("dry air", 295.0).status == "validated"
```

## Qualified intervals

A dispersion model's source supports it over one or more disjoint intervals of vacuum
wavelength, `intervals`; `wavelength_range` is their envelope. A wavelength outside
every interval warns (`ModelRangeWarning`); a NaN does not. Single-interval models use their stated range. A fit configured with OEP 3's `windows = true` has several,
the windows between the molecule's absorption bands, and between two of them the
warning names the gap's reason (`absorption`, `representation`, or `width`) (see
{doc}`fitting` and {doc}`OEP 3 <oeps/archive/oep-0003-infrared-windows>`).

```python
assert n2.intervals == (n2.wavelength_range,)
```

## Ownership

A model owns its storage: its arrays are copied at construction and made read-only, and
its mappings are copied into read-only mappings. Changing an array after passing it to
a model changes nothing, and writing to a model's array raises. A model can be pickled
and sent to another process, and keeps its digest and read-only storage there. The
digest is computed once and cached.

```python
assert not n2.covariance.flags.writeable
```

This guards against accidents; it is not a security boundary.

## Explicit parameters and overrides

A model can be built from a caller's parameters, checked exactly as a curated one is:
`IonisationRateModel.from_parameters`, `MolecularRateModel.from_parameters`,
`OscillatorModel.from_parameters` (poles in SI), `GasRamanResponse.from_lines`, and for
the solids `FormulaModel.sellmeier` and `TabulatedModel.from_arrays`. The curated
builders go through the same constructors, so the two paths cannot differ. The species
may be a plain name, which becomes the identity's material (OEP 1, section 7.4). A
custom species has no curated virial ratio or equation of state: its index takes the
ideal-gas density.

`model.with_parameters(**changes)` rebuilds a model with some of its **inputs**
changed, the names each class declares in `INPUTS`:

| Class | Inputs |
|---|---|
| `IonisationRateModel` | the channels (and `ionization_energy`, which moves them all), `orbital_l`, `structure`, `structure_coefficient`, `electrons`, `coulomb`, `static_correction`, `static_constants`, `channel_sum`, `averaging`, `uncertainty_factor` |
| `MolecularRateModel` | `form`, the energies, `structure_coefficients`, `charge`, `electrons`, `coulomb`, the correction, `vibrational_constant`, `averaging`, `channel_sum`, `uncertainty_factor` |
| gas dispersion models | `intervals`; an oscillator model's `strengths`, `energies`, `labels`, `offset_au`, `retardation_au`; a fit's `infrared` and `infrared_scale`; a formula's `terms`, a database entry's `coefficients` |
| solid models | `intervals`; a formula's `terms`; a table's `knots` and `index`; a fitted UDM's `values` |
| `GasRamanResponse` | `amplitude_scale`, `width_scale`, `width_law` |
| `SilicaRamanResponse` | `parameters` (the shape's) and `raman_fraction` |

A field derived from the inputs is refused by name (`DERIVED`): the molecular structure
factors and Hartree coefficient, a range that follows its intervals, a fit's
covariance, the Raman damping groups and density range, silica's `paired`. Everything
derived is recomputed by the constructor's own rules. The equation-of-state, transport,
air, hyperpolarizability and Kerr models have no inputs: a different model of those is
a curated choice, not an override.

An overridden model has the status `user`. Its identity records each change in
`modified`, and the curated model's digest in `base`. Its stated uncertainty does not
carry over: a rate's `uncertainty_factor` and a fit's covariance become NaN unless the
change gives them.

```python
from optmat.gas.ionisation_rates import ionisation_rate_model

rate = ionisation_rate_model("Ar", 800e-9)
changed = rate.with_parameters(electrons=1.5)
assert changed.identity.status == "user"
assert changed.identity.modified == (("electrons", "1.5"),)
assert changed.identity.base == rate.identity.digest
narrow = raman_response("N2", 295.0).with_parameters(width_scale=0.8)
assert narrow.identity.status == "user"
```

This page documents what [OEP 1](oeps/archive/oep-0001-propagation-interface.md)
proposed and implemented: its theory, decisions and phase reports are the record.
