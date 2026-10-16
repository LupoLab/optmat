# API reference

Use this reference for signatures, parameters, units, and return types. Start with
the {doc}`/quickstart` for working examples or the {doc}`/user-guide/index` for
the models and their validation.

## Core

The curated records, their vocabulary and materials, states and wavelengths,
and model identity. See {doc}`/user-guide/foundations` for conventions and
{doc}`/the-data/index` for the data records.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :recursive:

   optmat.conditions
   optmat.constants
   optmat.datasets
   optmat.exceptions
   optmat.identity
   optmat.materials
   optmat.vocabulary
```

## Gases

The user guide covers {doc}`gas dispersion and absorption </user-guide/linear>`,
{doc}`nonlinear and Raman response </user-guide/nonlinear>`, and
{doc}`ionisation </user-guide/ionisation>`.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :recursive:

   optmat.gas
```

## Solids

See {doc}`/user-guide/solids` for the common interface and the material-specific
models and their sources.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :recursive:

   optmat.solid
```
