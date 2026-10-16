# Introduction

optmat supplies optical properties for simulations of ultrafast and nonlinear optics.
It brings together curated measurements, published formulas, and fitted models for
gases, glasses, crystals, and metals, with their sources and limitations.

Start with the functionality overview to see which materials and responses are
available. Then install the package and work through the quick start. The examples
assume familiarity with Python and optics; they introduce optmat's interfaces and
conventions rather than the underlying subjects.

```{toctree}
:maxdepth: 1

overview
installation
../quickstart
curation
../about
../changelog
```

## Finding the detail

The {doc}`/user-guide/index` develops the conventions, models, and their validation
in reading order. Its chapters also stand as entry points for particular properties.
For example, follow the quick start's gas-index calculation with
{doc}`/models`, then {doc}`/gas_models` and {doc}`/accuracy` to understand the
recommended model and its uncertainty.

Use {doc}`/the-data/index` to inspect the records behind a result and the evidence
for curation decisions. The {doc}`/api/index` gives signatures, units, and return
types. To add data, modify a model, or work on the package itself, continue to the
{doc}`/developers/index`.
