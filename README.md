# optmat

[![Tests](https://github.com/LupoLab/optmat/actions/workflows/ci.yml/badge.svg)](https://github.com/LupoLab/optmat/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/LupoLab/optmat/graph/badge.svg)](https://codecov.io/gh/LupoLab/optmat)
[![Docs](https://readthedocs.org/projects/optmat/badge/?version=latest)](https://optmat.readthedocs.io)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![basedpyright](https://img.shields.io/badge/basedpyright-checked-42b983)](https://docs.basedpyright.com)
[![PyPI](https://img.shields.io/pypi/v/optmat)](https://pypi.org/project/optmat/)
[![Python](https://img.shields.io/pypi/pyversions/optmat)](https://pypi.org/project/optmat/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![REUSE](https://img.shields.io/badge/REUSE-compliant-green)](https://reuse.software)

**Curated, validated optical properties of materials, with the provenance of every
number.**

optmat provides the optical response of materials for simulations of ultrafast and
nonlinear optics, in gases, glasses, crystals and metals:
- the refractive index and its dispersion, from the extreme ultraviolet to the infrared;
- absorption, and the complex index across absorption thresholds and bands;
- the Kerr and third-harmonic coefficients;
- the Raman response of the molecular gases, dry air and fused silica;
- the strong-field ionisation rates of the noble gases and of molecules, and the
  parameters those rates are built from;
- two-photon continuum cross sections and curated bound-excitation records for noble gases.

Every value it gives is traceable to the paper, table and page it came from, and each
model identifies its sources, validity limits and available uncertainty evidence. It is the materials library
behind the pulse-propagation code [solprop](https://github.com/LupoLab/solprop), and does
not depend on it.

The manual is at [optmat.readthedocs.io](https://optmat.readthedocs.io), and the
repository is [LupoLab/optmat on GitHub](https://github.com/LupoLab/optmat).

## Why

Simulations of nonlinear pulse propagation are only as good as the material data they
are given. A formula or a single coefficient needs its measurement conditions,
conventions, range of validity and uncertainty to be used reliably.

Our motivation is to fix transcription errors, incorrect readings of measurement
conditions and inappropriate applications of models, and to evaluate data that have
not been systematically assessed. optmat provides reference implementations of published
fits and develops new fits from the best available data and analysis. The evidence,
choices and remaining uncertainties are available for inspection and the calculations
can be repeated.

[refractiveindex.info](https://refractiveindex.info) is an invaluable resource for the
optics community, and we do not seek to replace it. Our goals are different and
complementary: optmat concentrates on systematic evaluation of the underlying data,
explicit physical conventions, reproducible fits and validated implementations for
simulations.

## Data and modelling approach

- **Provenance of every number.** Values are taken exactly as printed: extracted by
  script from each paper's tables, or transcribed from its text with the page. The
  source record states:
  - the technique and the conditions of the measurement;
  - the reference state and density convention;
  - the stated uncertainty and its confidence level.

  A printed value is never edited. A misprint is corrected by an *erratum* that records
  the printed value, the correction, the reason and the evidence.
- **A public review register.** [`REGISTER.md`](REGISTER.md) records, entry by entry:
  - what was transcribed and how it was checked;
  - every erratum and every doubt;
  - every decision taken in fitting, and every question still open.

  Code, data and documentation cite it by entry (`TC-92`).
- **Physical models.** The dispersion of a gas is modelled as the dipole
  oscillator-strength distribution of its atoms or molecules:
  - resonance lines at their spectroscopic energies, and effective oscillators for the
    rest;
  - vibrational bands from HITRAN for the infrared-active molecules;
  - ab initio values and sum rules where they exist.

  The spectroscopic constraints support interpolation and extrapolation, but do not
  establish accuracy outside the model's qualified range. The manual states the
  evidence and limitations for each gas.

  For fused silica, the complex dielectric function is modelled from its excitations: interband
  and excitonic transitions, the Urbach tail and phonon bands, as in Franta et al.'s
  universal dispersion model. It is causal and passive by construction, with an exact
  derivative from the extreme ultraviolet to the far infrared. Other solids use
  published formulas or curated tables with explicit interpolation rules.
- **Density conventions.** Measurements are reduced with the density the source
  meant: real-gas (from equations of state) or ideal. They are converted through
  Lorentz–Lorenz with the refractivity virial coefficient and the magnetic term. At the
  level of parts per million these corrections bring precision metrology and ab initio
  theory into agreement.
- **Uncertainty and sensitivity.** Fits store their covariance. Gas dispersion and
  hyperpolarizability fits are refitted
  without each of its datasets, to show which data it really rests on. Measurements not
  used in a fit are kept as validators and compared with it. The documentation's tables
  of accuracy are checked against the stored fits by the test suite.
- **Process-aware nonlinearity.** The second hyperpolarizability is modelled as a
  function of all four frequencies, so one model gives the Kerr effect, third-harmonic
  generation and four-wave-mixing channels within the model's stated support. An exact
  time-domain form of the model is derived and tested.
- **Responses checked against experiment and exact theory.** The noble gases'
  ionisation rates are validated against numerical solutions of the Schrödinger
  equation and against measured yields. The Raman lines are checked against measured
  gains.

## What is there now

- **761 source records:** 8467 datasets and 197,754 values, 259 errata and 780
  published formulas.
- **15 gases:** He, Ne, Ar, Kr, Xe, H2, D2, N2, O2, CO, CO2, CH4, SF6, N2O and O3.
  - **A fitted model of the refractive index for every gas:**
    - typically from the vacuum ultraviolet to the mid-infrared;
    - for helium and argon, at the accuracy of ab initio theory (0.1 ppm and 3 ppm at
      800 nm);
    - for ozone, a complex index built from its measured absorption spectrum.
  - **The extreme ultraviolet:** the complex index of the noble gases, across their
    absorption thresholds.
  - **The third-order response:** a model of the second hyperpolarizability for every
    gas, with the conversions from every convention in the literature.
  - **The Raman response:**
    - the Raman lines of N2, O2, H2, D2, CO2, CH4 and SF6, and N2O's vibrational
      bands, with per-line widths as functions of density and composition;
    - the response functions of the gases and of dry air, checked against measured
      Raman gains.
  - **Ionisation:**
    - the parameters of every rate in use, from the primary data;
    - validated rates of the noble gases, with their accuracy stated;
    - a rate for the molecules, set by all-electron calculations of H2 and N2, with
      the forms of the literature beside it.
- **Solids:** 95 of the records hold data of glasses, crystals and metals: 2027
  datasets for 163 materials:
  - the linear index and the band gap;
  - n2 and the Raman response;
  - two-photon and multiphoton absorption;
  - ionisation and laser damage.
- **Fused silica, the first solid modelled:**
  - its complex index from 25 nm to 125 µm: Franta et al.'s universal dispersion model
    refitted to the absolute refractometry, which it follows to about 2e-5 in the
    transparent window;
  - its n2 from 350 nm to 3.5 µm, with the Raman fraction set by the measured Raman
    gain and the dependence on pulse duration;
  - its Raman response function.
- **Crystalline silicon:** its complex index from 31 nm to 143 µm at 20 °C, Franta
  et al.'s universal dispersion model refitted to the absolute refractometry (to
  1.4e-4 from 1.2 to 5.5 µm) and the measured absorption edge.
- **The window crystals and N-BK7:** CaF2, MgF2 and sapphire (both rays), BaF2 and
  SCHOTT's N-BK7, each its source's Sellmeier formula made absolute (the prism and
  catalogue indices are relative to air or nitrogen, which the formulas in common use
  ignore: about 4e-4 in n).
- **The metals:** silver, gold and copper from the soft x-ray range (3-4 nm) to
  12.4 µm (copper to 56 µm), aluminium from 0.8 nm to 200 µm, as curated tables joined
  by range.
- **Every published formula** of the curated sources, available by name for
  comparison.

The [roadmap](https://optmat.readthedocs.io/en/latest/roadmap.html) describes what comes
next:
- the models of the other solids, and their two-photon absorption, ionisation and
  damage;
- the molecules in the extreme ultraviolet;
- line mixing in the dense molecular gases, and the remaining reference rates.

## Installation

optmat needs Python 3.12 or later:

```sh
pip install "optmat[gas]"
```

Or, with [uv](https://docs.astral.sh/uv/), from a clone of the repository:

```sh
git clone https://github.com/LupoLab/optmat.git
cd optmat
uv sync --locked --no-dev --extra gas
uv run --no-sync python
```

The `gas` extra supplies CoolProp for real-gas densities. The solids
need no extra. For pip and platform-specific instructions, see
[Installation](https://optmat.readthedocs.io/en/latest/introduction/installation.html).
Use `uv sync --locked` for the development environment.

## A first look

```python
from optmat.conditions import State
from optmat.gas.gases import Gas
from optmat.gas.models.registry import default_model, refractivity
from optmat.gas.nonlinear import hyperpolarizability_models

# n - 1 of helium at 800 nm, 30 bar and 20 °C, from the recommended model
refractivity("He", 800e-9, State(30e5, 293.15))

# the model behind it, with its stated range and citation
default_model("He")

# argon's Kerr coefficient n2 at 800 nm and 1 bar, and its third-harmonic gamma (au)
argon = hyperpolarizability_models()[Gas.Ar]
argon.n2(800e-9, State(1e5, 293.15))
argon.gamma_au("thg", 800e-9)
```

The Raman response and the ionisation rate:

```python
from optmat.gas.ionisation_rates import ionisation_rate_model
from optmat.gas.linewidths import AMAGAT
from optmat.gas.raman_response import raman_response

# the Raman response h(t) of hydrogen at 295 K and 2 amagat, 1 ps after the pulse
hydrogen = raman_response("H2", 295.0)
hydrogen.response(1e-12, 2 * AMAGAT)

# argon's ionisation rate at 800 nm, averaged over a cycle of a 1e14 W/cm^2 field (1/s)
ionisation_rate_model("Ar", 800e-9).cycle_averaged_rate(2.74e10)
```

Fused silica:

```python
from optmat.solid.kerr import kerr_models
from optmat.solid.models.registry import default_model

silica = default_model("SiO2")  # fit/SiO2, at 20 degrees C
silica.refractive_index(120e-9)  # complex n + ik in the vacuum ultraviolet
silica.gvd(800e-9)  # group-velocity dispersion, s^2/m

# the long-pulse n2 at 1030 nm, and what a 50 fs pulse measures (m^2/W)
kerr = kerr_models()["kerr/SiO2"]
kerr.n2(1030e-9)
kerr.n2_effective(1030e-9, 50e-15)
```

The manual has five entrances:

- [Introduction](https://optmat.readthedocs.io/en/latest/introduction/index.html): capabilities, installation and the quick start.
- [User guide](https://optmat.readthedocs.io/en/latest/user-guide/index.html): conventions, models, interpretation, accuracy and limitations.
- [The data](https://optmat.readthedocs.io/en/latest/the-data/index.html): the generated source catalogue, record format, review register and gallery.
- [Developers guide](https://optmat.readthedocs.io/en/latest/developers/index.html): contribution, extraction, fitting and verification workflows.
- [API](https://optmat.readthedocs.io/en/latest/api/index.html): signatures, units and return values.


## Contributing

Corrections, new data and new materials are welcome. Read
[CONTRIBUTING.md](CONTRIBUTING.md) first: it explains how data enter the package and
what a contribution must record to keep every value traceable. Major changes start as
enhancement proposals (OEPs, `docs/oeps/`), which are reviewed before they are
implemented and keep the record of how each feature was built and verified.

## AI-assisted development and verification

AI agents assist with literature curation, implementation, tests and documentation
under John C. Travers's direction. He decides which data, models and conventions are
accepted. Source-specific extraction scripts, checks against the primary papers,
numerical references and executable examples support that work; the review register
records corrections, decisions and unresolved questions.

Read [About optmat](https://optmat.readthedocs.io/en/latest/about.html) for responsibility
and verification, and the [Developers guide](https://optmat.readthedocs.io/en/latest/developers/index.html)
for the technical procedures.

## Licence

The code, the source records and the fits are [MIT](LICENSE); the manual's prose and
figures are [CC-BY-4.0](LICENSES/CC-BY-4.0.txt). optmat is self-contained: it ships the
published data its models use, transcribed into `src/optmat/data/values/` or copied
unchanged into `src/optmat/data/external/`. Those data belong to their authors and
publishers, are cited in every case, and are included for scientific research under
the [third-party data notice](LICENSES/LicenseRef-optmat-third-party.txt), which also
says how a rights holder can ask for their terms to be applied. Please cite the
original sources, as well as optmat, when you use the values.
[`REUSE.toml`](REUSE.toml) gives the copyright and licence of every file.
