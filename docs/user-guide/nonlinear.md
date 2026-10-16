# Nonlinear and Raman response

The third-order response depends on the optical process and its frequencies. Raman
response also requires a distinction between electronic, vibrational, and rotational
parts, and between a line's strength and the response in time.

The {doc}`/quickstart` gives examples of a Kerr coefficient and a Raman response.
Continue below for their conventions, physical models, support, and validation.
Use {doc}`/conventions` for shared units and states. The next chapter is
{doc}`ionisation`.

```{toctree}
:maxdepth: 1

../nonlinear
nonlinear-conventions
nonlinear-models
nonlinear-support
../nonlinear_accuracy
nonlinear-extrapolation
../raman
../raman_response
raman-validation
```

The {doc}`nonlinear gallery </auto_gallery/nonlinear/index>` illustrates these
responses. Reference interfaces include {mod}`optmat.gas.nonlinear`,
{mod}`optmat.gas.raman`, and {mod}`optmat.gas.raman_response`.
