# Developers guide

This guide is for work on optmat: its curated records, models, implementation, and
documentation. To use the existing package in a calculation, start with
{doc}`/introduction/installation` and the {doc}`/user-guide/index`.

## Development environment

Clone the [repository](https://github.com/LupoLab/optmat) and install its locked
development environment with [uv](https://docs.astral.sh/uv/):

```sh
git clone https://github.com/LupoLab/optmat.git
cd optmat
uv sync --locked
uv run poe check
```

The development environment includes the gas extra, testing tools, and documentation
dependencies. The full gate checks style, types, docstrings, spelling, tests with
statement and branch coverage, and the manual including its executable gallery.

Read the [contributor instructions](https://github.com/LupoLab/optmat/blob/main/CONTRIBUTING.md)
and [working rules](https://github.com/LupoLab/optmat/blob/main/AGENTS.md) before
making a change. They describe the additional checks required for changes to data,
fits, packaging, and public interfaces.

## Choosing the next step

For a data contribution, read {doc}`/data` and {doc}`/register`, then follow
{doc}`/contributing`. For model work, first read the relevant scientific chapter
in the {doc}`/user-guide/index`; {doc}`/fitting` explains the gas dispersion fits.
A major change follows the enhancement-proposal process. The roadmap records
accepted proposals, drafts, and remaining work.

```{toctree}
:maxdepth: 1

testing
../contributing
data
fitting
documentation
register
../oeps/index
../roadmap
```

## Package architecture

optmat has a medium-independent core and subpackages for gases and solids:

- **Core:** curated records, quantity vocabulary, materials, conditions, physical
  constants, model identity, and exceptions.
- **Gases:** densities, refractive-index models and fitting, EUV response,
  third-order response, Raman lines and response functions, ionisation, mixtures,
  and transport.
- **Solids:** complex-index models and tables, fitting, and fused silica's Kerr and
  Raman response.

The core does not import a medium's subpackage. The {doc}`/api/index` lists the
modules, and {doc}`/propagation` describes the interfaces used by propagation codes.
