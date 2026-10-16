# Quick start

Follow {doc}`introduction/installation` to install optmat with its `gas` extra.
The examples below introduce the main interfaces; run the Python blocks in order
in one session. Every block runs as part of the test suite and the installation
checks.

Wavelengths are vacuum wavelengths in metres, pressures are in pascals, and
temperatures are in kelvin. The {doc}`introduction/overview` describes the available
properties, and the {doc}`user-guide/index` develops the models and their validation.

## The refractive index of a gas

Every gas has a recommended model, chosen for stated reasons. The one-call entry point
gives $n - 1$ at any wavelength and state. Wavelengths are vacuum wavelengths in metres;
states are a pressure in pascals and a temperature in kelvin.

```python
from optmat.conditions import State
from optmat.gas.models.registry import refractivity

state = State(1e5, 293.15)  # 1 bar, 20 °C
refractivity("He", 800e-9, state)  # n - 1 of helium at 800 nm
refractivity("Ar", 800e-9, State(30e5, 293.15))  # argon at 30 bar
```

The density at the state is the real-gas density, and $n$ follows from Lorentz–Lorenz
with the refractivity virial term. The page {doc}`conventions` explains both.

## Choosing and inspecting a model

Every gas has several models: the fit made here (`fit/<gas>`), each published formula of
the curated sources (`<source>/<formula>`), and the entries of the refractiveindex.info
database, for comparison.

```python
import numpy as np
from optmat.gas.models.registry import default_model, model, models

sorted(models("He"))  # every model of helium, by name
helium = default_model("He")  # the recommended one
helium.name, helium.wavelength_range, helium.citation

wavelengths = np.array([500e-9, 800e-9])
helium.refractivity(wavelengths, state)
model("borzsonyi_2008/sellmeier_He").refractivity(wavelengths, state)
```

A model warns (`ModelRangeWarning`) when it is evaluated outside the range its data
support: Börzsönyi et al.'s helium formula at 1.6 µm, for instance, which they fitted
from 400 to 1000 nm.

Continue with {doc}`models` for model selection, {doc}`gas_models` for the
recommended models, and {doc}`identity` for status, qualified intervals, and the
information to save with a result.

## How well a model is known

A fitted model carries its covariance. The relative one-sigma uncertainty of the
group-velocity dispersion ($k = 2$) or third-order dispersion ($k = 3$) is independent
of the density in the dilute limit:

```python
from optmat.gas.fitting import HARTREE_NM
from optmat.gas.models.fitted import fitted_models

argon = fitted_models()["fit/Ar"]
photon_energy = np.array([HARTREE_NM / 800.0])  # 800 nm, in hartree
argon.dispersion_sigma_relative(photon_energy, 2)
```

The page {doc}`accuracy` gives these uncertainties for every gas, with the leave-one-out
spread and a working estimate.

## Absorbing gases and the extreme ultraviolet

Ozone's model is complex: its measured absorption enters by Kramers–Kronig. The noble
gases have a complex model across their absorption thresholds.

```python
from optmat.gas.euv import euv_model
from optmat.gas.models.registry import refractive_index

refractive_index("O3", 600e-9, state)  # complex n of ozone in the Chappuis band
euv_model("Ar").refractive_index(60e-9, state)  # argon at 60 nm, above threshold
```

See {doc}`euv` for the noble-gas models and their sources, and {doc}`gas_models`
for ozone and the other recommended dispersion models.

## The third-order response

The second hyperpolarizability of each gas is a model of the process and the
frequencies. It gives the Kerr coefficient and the third-harmonic coefficient;
{doc}`nonlinear` explains where the fitted response is supported:

```python
from optmat.gas.gases import Gas
from optmat.gas.nonlinear import hyperpolarizability_models

argon_gamma = hyperpolarizability_models()[Gas.Ar]
argon_gamma.gamma_au("kerr", 800e-9)  # gamma for SPM and XPM, atomic units
argon_gamma.gamma_au("thg", 800e-9)  # gamma for the third harmonic of 800 nm
argon_gamma.n2(800e-9, state)  # n2 at 1 bar, m^2/W
```

The page {doc}`nonlinear` explains the conventions and the model.

## The Raman response

Evaluate the molecular gases' Raman lines and the response function built from
them at a specified temperature and density:

```python
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman import raman_lines
from optmat.gas.raman_response import raman_response

lines = raman_lines("H2", 295.0, 2 * AMAGAT)  # at 295 K and 2 amagat
lines.branch, lines.j, lines.frequency, lines.fwhm  # each line, and its width (Hz)

air = raman_response("dry air", 295.0)
air.response(np.array([0.5e-12, 1e-12]), AMAGAT)  # h(t), m^5 V^-2 s^-1
```

The pages {doc}`raman` and {doc}`raman_response` give the convention, the sources and
the checks against measured gains.

## Ionisation

The parameters every rate in use is built from, and optmat's validated rates:

```python
from optmat.gas.ionisation import ionisation_parameters
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.molecular_ionisation_rates import molecular_rate_model

ionisation_parameters("Kr").ionization_energy  # J, NIST ASD
xenon = ionisation_rate_model("Xe", 800e-9)
xenon.cycle_averaged_rate(2.74e10)  # 1/s, at the amplitude of 1e14 W/cm^2
xenon.uncertainty_factor  # the factor within which it is expected to be right
molecular_rate_model("N2", 800e-9).instantaneous_rate(2e10)  # along a field, V/m
```

The pages {doc}`ionisation`, {doc}`ionisation_rates` and
{doc}`molecular_ionisation_rates` give the forms and their validation.

## Solids

Fused silica, the window crystals, N-BK7 and the metals have recommended models. A
solid's index is at its reference temperature, and needs no state:

```python
from optmat.solid.kerr import kerr_models
from optmat.solid.models.registry import default_model

silica = default_model("SiO2")
silica.refractive_index(800e-9), silica.gvd(800e-9)  # complex n, and beta_2 in s^2/m
kerr_models()["kerr/SiO2"].n2(800e-9)  # m^2/W
default_model("MgF2", "o").refractive_index(1030e-9)  # a crystal needs its axis
default_model("Au").refractive_index(800e-9)  # a metal: n + ik
```

The pages {doc}`silica`, {doc}`silicon`, {doc}`crystals` and {doc}`metals` give the models and what
each rests on.

## The data behind the models

Every curated value is available with its provenance:

```python
from optmat.datasets import datasets, sources
from optmat.gas import reduction

refractivities = datasets("Ar", quantity="refractivity")
dataset = refractivities[0]
dataset.key, dataset.state, dataset.scaled_values()  # as printed, times the scale
sources()[dataset.source].citation  # the paper
reduction.polarizability(dataset)  # the polarizability it implies
```

The page {doc}`data` describes the records, and {doc}`register` how each was checked.
For an overview of the route from a paper to a model, read
{doc}`introduction/curation`.

## Continuing with the user guide

Read {doc}`user-guide/foundations` next for the shared conventions and model
identity. For a particular response, use the chapters in {doc}`user-guide/index`.
If you are integrating optmat into a simulation, {doc}`propagation` brings the
interfaces together; {doc}`api/index` supplies the detailed reference.
