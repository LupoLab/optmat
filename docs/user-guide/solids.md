# Solids

Solid models provide the complex index at their stated reference conditions. Read
the common interface and conventions first, then the chapter for the material of
interest. For an anisotropic crystal, the optical axis is part of the model choice.

The {doc}`/quickstart` evaluates fused silica, a crystal, and a metal. These chapters
explain their sources, ranges, and limitations, including silica's nonlinear
response. The next step in the guide is {doc}`/propagation`.

```{toctree}
:maxdepth: 1

../solid
../silica
../silicon
../crystals
../metals
```

See the {doc}`solid gallery </auto_gallery/solids/index>` and the reference
interfaces {mod}`optmat.solid.models.registry`, {mod}`optmat.solid.kerr`, and
{mod}`optmat.solid.bundle`.
