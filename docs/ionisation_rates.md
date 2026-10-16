# Ionisation rates

The rate at which an intense, linearly polarised field ionises a noble-gas atom:
{mod}`optmat.gas.ionisation_rates`. It is built from published theory and published
constants, and validated against all-electron calculations, time-dependent Schrödinger
calculations and measured yields ({mod}`optmat.gas.ionisation_validation`). The
parameters it rests on are those of {doc}`ionisation`.

```python
import numpy as np

from optmat.gas.ionisation_rates import ionisation_rate_model

field = 2.74e10 * np.cos(np.linspace(0.0, 2 * np.pi, 101))  # one cycle, V/m

argon = ionisation_rate_model("Ar", 800e-9)
argon.cycle_averaged_rate(2.74e10)  # 1/s, at the amplitude of 1e14 W/cm^2
argon.instantaneous_rate(field)  # 1/s, along a carrier-resolved field in V/m
argon.keldysh(2.74e10), argon.barrier_suppression_field, argon.maximum_field
argon.electrons, argon.channels, argon.static_constants  # what it is made of
argon.citation, argon.sources, argon.reasons  # and where each part comes from
```

The gases are He, Ne, Ar, Kr and Xe. The field is in V/m and the rate in s⁻¹; the
amplitude of an intensity $I$ is $F = (2I/\varepsilon_0 c)^{1/2}$.

Read the {doc}`/quickstart` for the shorter first calculation. Continue with
{doc}`user-guide/ionisation-theory` for the physical definitions, or
{doc}`user-guide/ionisation-accuracy` for the supported range and uncertainty.

## The models

<!-- ionisation-rate-models:start -->
| Model | Coulomb factor | Static correction | Coefficient | Fine structure | N_eff |
|---|---|---|---|---|---|
| `default` | frequency | tong-lin | tong_2002 | yes | 1, 1.35, 1.47, 1.5 |
| `ppt_larochelle` | static | none | hartree | no | 1 |
| `ppt_ilkov` | ilkov | none | hartree | no | 1 |
| `ppt_popruzhenko` | popruzhenko | none | tong_2002 | no | 1 |
| `one_electron` | frequency | tong-lin | tong_2002 | no | 1 |
| `default_static` | static | tong-lin | tong_2002 | yes | 1, 1.47 |
| `one_electron_static` | static | tong-lin | tong_2002 | no | 1 |
<!-- ionisation-rate-models:end -->

`default` is the model above. `one_electron` is the default without $N_\text{eff}$
and the fine structure, which is what a one-electron calculation in Tong & Lin's
potentials computes. `default_static` and `one_electron_static` are the same two with
the static Coulomb factor, the default until TC-848; the comparisons of the register
before then were made with them. The `ppt_*` models are the rates of **one electron**
as their sources give them, for comparison.

## Validation

The validation separates static normalisation, frequency dependence, and the
atom's effective electron count. Some references set parameters: helium's electron
count uses Parker et al., the frequency-dependent Coulomb correction is fitted to
one-electron argon calculations, and Kr and Xe electron counts use Wallace et al.
The chapters identify those choices alongside independent tests.

Tables are generated from {mod}`optmat.gas.ionisation_validation` and checked against
the stored results. `default` and `one_electron` use the frequency-dependent Coulomb
factor adopted in TC-848; `default_static` and `one_electron_static` preserve the
older static-factor calculation. Historical comparisons are identified explicitly.

<span id="the-model"></span>
<span id="three-parts-each-tested-against-what-isolates-it"></span>
<span id="two-rates"></span>
<span id="the-energy-per-event"></span>
<span id="what-is-left-out-and-why"></span>

See {doc}`/user-guide/ionisation-theory`.

<span id="the-effective-number-of-electrons"></span>
<span id="static-rates"></span>

See {doc}`/user-guide/ionisation-static`.

<span id="helium-in-a-laser-field"></span>
<span id="all-electron-rates-at-390-nm"></span>

See {doc}`/user-guide/ionisation-laser`.

<span id="one-electron-calculations-of-argon-800-nm-6-fs-and-wallace-et-al-s-5-5-fs-400-800-and-1030-nm-10-fs-400-and-1030-nm-30-fs"></span>
<span id="one-electron-calculations-of-argon-800-nm-6-fs-and-wallace-et-als-55-fs-400-800-and-1030-nm-10-fs-400-and-1030-nm-30-fs"></span>
<span id="one-electron-calculation-of-argon-800-nm-42-fs"></span>
<span id="one-electron-calculations-of-hydrogen"></span>
<span id="the-frequency-dependent-coulomb-factor"></span>

See {doc}`/user-guide/ionisation-one-electron`.

<span id="measured-yields"></span>
<span id="multiphoton-saturation-at-248-6-nm"></span>
<span id="circular-polarisation"></span>

See {doc}`/user-guide/ionisation-measured`.

<span id="accuracy"></span>
<span id="open-questions"></span>
<span id="luna-jl-and-solprop"></span>

See {doc}`/user-guide/ionisation-accuracy`.
