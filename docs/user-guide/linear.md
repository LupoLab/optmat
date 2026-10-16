# Gas dispersion and absorption

Start with {doc}`/models` to select and evaluate a gas model. The
{doc}`/quickstart` contains working examples of the refractive index, inspecting a
model, and evaluating an absorbing gas. Read {doc}`/conventions` for the relation
between polarizability, density, and index.

The chapters below explain the fitting method, what each recommended model rests
on, and its accuracy. They then extend the discussion to the extreme ultraviolet,
mixtures, and transport properties. In the guide's reading order, {doc}`nonlinear`
follows these linear properties.

```{toctree}
:maxdepth: 1

../models
../fitting
../gas_models
../accuracy
../euv
../mixtures
../transport
```

Inspect the {doc}`dispersion gallery </auto_gallery/dispersion/index>` and
{doc}`accuracy gallery </auto_gallery/accuracy/index>` alongside the chapters.
The main reference interfaces are {mod}`optmat.gas.models.registry`,
{mod}`optmat.gas.euv`, {mod}`optmat.gas.mixture`, and {mod}`optmat.gas.transport`.
