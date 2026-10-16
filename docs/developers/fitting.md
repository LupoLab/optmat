# Running fits and robustness checks

For the scientific method and interpretation, read {doc}`/fitting` and
{doc}`/nonlinear`. This page describes the maintenance workflow. A change to a
fit, its code or its input data requires the checks below and register entries for
selection, priors, scales and ranges, with the numbers that motivated each decision.

## Dispersion fits

Configurations are in `src/optmat/data/gas/fits/`; fitted records are in
`src/optmat/data/gas/models/fitted/`.

```sh
uv run python tools/fit/run.py <gas>
uv run python tools/fit/robustness.py <gas>
```

The first command stores the fitted parameters, covariance and validation. The second
refits without each dataset, starting from the full solution, and stores the
leave-one-out record under `data/gas/models/robustness/`.

Phase-matching observations hold their partner at its recommended model. Refit the
affected partners first: helium and argon; then neon, nitrogen, hydrogen, krypton and
xenon; then oxygen and deuterium. A changed partner also requires refitting the gases
that depend on it. The configurations identify which ratios are fitted or validated.

## Infrared configuration

For an infrared-active molecule, `[model] infrared_max_wavenumber` adds the gas's
HITRAN bands up to that wavenumber as fixed terms (not fitted; see the data page and
`optmat.gas.infrared`): measured refractivities and dielectric constants include them, while
electronic theory is marked `infrared = false` and Cauchy moments are always electronic.
The fixed terms are stored after the fitted ones and are outside the covariance.

Two keys read the gas's derived tables instead (OEP 3; `data/external/hitran/<gas>_*.csv`,
{doc}`/data`). The CO, CO2, CH4 and N2O fits use these band rules and windows
(TC-797):

- `[model] band_nodes = N` takes every band system's rule, at most `N` nodes each, as
  the infrared terms, in place of `infrared_max_wavenumber`'s band poles (the two are
  exclusive, and `N` must reach the tables' own node count). The rotational system of a
  polar molecule becomes a term class of its own: it enters refractive observations
  with the infrared terms and is left out of electric ones (static dielectric data)
  and of the moments.
- `[model] windows = true` makes the gas's derived windows the model's qualified
  intervals, with each gap's reason.

A fit that uses either records `intervals_nm`, `gap_reasons`, `n_rotational`,
`band_nodes` and the SHA-256 of the tables it read (`derived_sha256`); a fit that uses
neither records nothing new. Such a fit gives `representation_bound_au(omega, order)`,
the bound on its infrared terms' representation error for order 0, 2 or 3, in every
qualified interval, and NaN outside them.

## Compact models

`data/gas/compact/<gas>.toml` holds a reduction's settings: the parent,
the number of poles per absorbing region (in ascending energy; a region given as many
poles as it holds keeps the fit's own), an optional `longest_nm` where the compact
intervals end, the grid, the tolerance and the random starts with their seed.
`uv run python tools/fit/compact.py <gas>` derives the poles from the stored fit, prints
the residuals and the reproduction of every dataset and writes
`data/gas/models/compact/<gas>.json`; `--fewer` also shows that one pole fewer fails.
The stored record carries the SHA-256 of the stored fit and of both configurations, so a
default-lane test fails when the fit changes without its compact model; the `fits` lane
re-derives every compact model. Run it after every refit of CO, CO2, CH4, N2O or SF6.


## Third-order fits

### Configuration

Each gas has a configuration, `data/gas/nonlinear/<gas>.toml`, read by
{func}`~optmat.gas.nonlinear.load_nonlinear_config`:

| Key | Meaning |
|---|---|
| `gas`, `description` | the species and what the fit is for |
| `range_nm` | the fundamental vacuum-wavelength range over which the model is declared valid |
| `b_prime_prior` | `[value, sigma]` prior on $B'$; absent fixes $B' = 0$ |
| `a_prior`, `b_prior` | optional priors on $A$ and $B$: a dispersion borrowed from a similar molecule (N2O) |
| `[[fitted]]` | the datasets fitted, each a key or a table of options |
| `validators` | the datasets compared with the fit, not fitted |

The options of a fitted dataset ({class}`~optmat.gas.nonlinear.NonlinearOptions`):

| Option | Meaning |
|---|---|
| `key` | `source/dataset` |
| `extra_relative` | a relative uncertainty added in quadrature to every value: the allowance for a calculation, or for a value without a stated uncertainty (every fitted value needs one) |
| `free_scale` | fit the dataset's level without a prior: only its shape informs the fit |
| `systematic` | fit the level with this relative prior (one sigma) |
| `scale_group` | datasets naming the same group share one scale: a calculation's processes, whose relative dispersion informs $B'$ whatever its level |
| `range_nm` | the fundamental-wavelength window used (a calculated THG, for example, only above 650 nm) |
| `rows` | only these rows (one method of a table of methods) |
| `subtract_dataset`, `subtract_formula` | the $\gamma^{vR}$ subtracted to leave the electronic part |
| `factor` | a unit factor the conversion does not know, with its reason stated in the file |

### Python fitting interface

To refit after changing a configuration:

```python
from optmat.gas.nonlinear import (
    fit_hyperpolarizability,
    load_nonlinear_config,
    validation,
)

result = fit_hyperpolarizability(load_nonlinear_config("Ar"))
result.parameters, result.chi2, result.dof, result.scales
validation(result)
```

### Running and storing

```text
uv run python tools/fit/nonlinear.py <gas> [<gas> ...]
uv run python tools/fit/nonlinear_robustness.py <gas> [<gas> ...]
```

- **The first** fits the gas, prints the parameters, the per-dataset $\chi^2$ and scales
  and the validation, and writes `data/gas/models/nonlinear/<gas>.json`:
  - the parameters and their covariance;
  - the largest fitted $W_2$;
  - the scales, $\chi^2$ and degrees of freedom;
  - the validation;
  - the SHA-256 of the configuration.

  A test re-runs every stored fit and compares.
- **The second** refits without each fitted dataset or scale group, and writes
  `data/gas/models/nonlinear_robustness/<gas>.json` (the leave-one-out records of the
  accuracy page).

After a refit, re-derive the gas's support record ({doc}`/user-guide/nonlinear-support`):

```text
uv run python tools/fit/nonlinear_support_record.py <gas> [<gas> ...]
```

## Verification and documentation

After all affected fits and derived records have been regenerated, run:

```sh
uv run poe test-fits
uv run poe check
```

Update the numerical documentation tables and their narrative when a finding changes.
The {doc}`documentation` page describes table generation and freshness tests. A
hyperpolarizability refit also requires its support record, even when the best-fit
parameters change only slightly. Dispersion refits of CO, CO2, CH4, N2O or SF6 also
require the compact-model derivation above. These derived records include hashes so
that stale results cannot silently accompany a changed fit.
