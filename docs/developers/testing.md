# Development checks

Run the full gate before considering a change complete:

```sh
uv run poe check
```

It runs formatting and lint checks, basedpyright, numpydoc, spelling, the copyright and
licence check (`reuse lint`), the tests under the coverage gate, and Sphinx with
warnings treated as errors. The gallery executes during the build. CI
(`.github/workflows/ci.yml`) runs the same gate on Linux and macOS for every push and
pull request, the tests on Windows, and the tests on Python 3.12 with the lowest
versions of the package's direct dependencies, the floors in `pyproject.toml`.
The commands and requirements in
[AGENTS.md](https://github.com/LupoLab/optmat/blob/main/AGENTS.md) are the project policy.

## Additional verification

| Change | Required check |
|---|---|
| Extracted data | `PAPER_LIBRARY=<path> uv run poe test-library` |
| Fits, fitting code or data used by a fit | Refit, regenerate robustness and derived records as described in {doc}`fitting`, then `uv run poe test-fits` |
| Packaging, data layout or public interfaces | `uv run poe test-release` |
| Documentation | `uv run poe docs`, relevant snippet and generated-table tests, and the full gate |

A proposal-only change uses the shorter gate specified in the working rules
(`AGENTS.md`): diff whitespace, spelling,
proposal and register tests, and the documentation build. This exception does not
apply when implementation changes accompany a proposal.

Tests run in parallel in the gate and must be independent of order and share no
writable output. Use `tmp_path` for test output. A focused invocation such as
`uv run pytest tests/test_docs_snippets.py` runs serially.

## Numerical regression checks

`tests/test_golden_models.py` compares sampled model outputs with
`tests/data/golden_models.json` at a relative tolerance of `1e-12`, with no absolute
floor. This allows arithmetic rounding, including for very small SI quantities.
Versions, classes and model metadata are also checked through their identities.

Gas Raman responses derive width coefficients by solving a small linear system.
Different floating-point arithmetic can change their exact identity digest while
preserving outputs within this tolerance. Their canonical identity records are stored
once per response in `tests/data/golden_raman_identities.json`. The test checks each
record's exact SHA-256, compares floating fields at `1e-12` relative tolerance, and
requires exact versions, classes, line identifiers, damping-group labels, source
references and other discrete metadata. Other model families retain exact digest
comparisons. The public {doc}`/identity` contract remains unchanged.

Regenerate both golden files after a reviewed model change:

```sh
uv run python tools/golden_models.py
```

Review the numerical and identity changes together, update the changelog and raise
the relevant `EVALUATION_VERSION` when evaluation changes. A failed digest check
requires diagnosis: it does not by itself establish a scientific change.

`tests/test_nonlinear_fits.py` compares each fresh hyperpolarizability fit with its
stored model using complementary numerical checks:

- The configuration digest, dataset selection, point counts and degrees of freedom
  agree exactly.
- Parameter displacements are at most `1e-4` of each stored standard uncertainty;
  fixed parameters agree exactly. A near-zero coefficient is assessed on its
  uncertainty scale, rather than relative to its small fitted value.
- Predictions agree to `1e-7` relative at every fitted observation and on a grid over
  the configured wavelength range, including its endpoints, for static, dc Kerr,
  Kerr, ESHG and THG response. Dataset scales agree to the same relative tolerance.
- The full least-squares objective, including parameter priors and each shared scale
  prior once, differs by at most `1e-8 * max(1, recorded objective)`. Stored and fresh
  data chi-squared values must also reproduce their respective parameters.

The fitter uses default stopping tolerances of `1e-8`. Near a least-squares minimum,
parameter resolution scales as the square root of objective resolution in uncertainty
units; this motivates the `1e-4` parameter bound. Convergence experiments with default
finite differences, three-point differences and an analytic Jacobian gave parameter
displacements below `7e-6` standard uncertainties and relative prediction changes below
`1.4e-8` at the fitted points. The prediction bound allows modest arithmetic variation
beyond those experiments while remaining much smaller than model uncertainty. These
are regression tolerances, not scientific error bars or guarantees on every platform.

Mutation tests verify that these checks reject changed model parameters, versions,
metadata, outputs, dataset selections, scales and objective values, while accepting
round-off. The stored models and runtime identity hashes are not rounded or rewritten
to make the tests pass. Broader platform verification remains in {doc}`/roadmap`.

`tests/solid/test_solid_fitting.py` compares each fresh solid fit (`fit/SiO2`,
`fit/Si`) with its stored record in the same way (TC-1089):

- The configuration digest, degrees of freedom, offsets and per-dataset point counts
  agree exactly.
- The free parameters, group scales and offsets move by at most `1e-2` of their stored
  standard uncertainties.
- On a grid over the model's range, n agrees to `1e-5` relative and k to `1e-4`
  relative wherever k exceeds `1e-12`.
- chi-squared differs by at most `1e-8 * max(1, recorded chi-squared)`.

These fits stop on the optimiser's `1e-8` tolerances at a point that depends on its
path along flat directions of the objective. Three settings of the unknowns' scaling
moved silica's parameters by up to `3.3e-3` standard uncertainties, n by `1.9e-6`
relative (in the extreme ultraviolet; `2e-9` in the transparent window) and chi-squared
by `2e-9` relative. The bounds are three to five times these. Mutation tests check that
a 0.05-sigma parameter shift, a changed objective, degrees of freedom or dataset
selection, and a changed group scale are rejected.

`tests/test_stored_fits.py` compares each fresh gas dispersion fit with its stored
record on its predictions, not its parameters (TC-1090). The effective oscillators
are nearly degenerate, and refits move single parameters by many of their stored
standard uncertainties while the model hardly moves:

- The degrees of freedom and per-dataset point counts agree exactly.
- Within each qualified interval, alpha agrees to `0.05` of the stored model's own
  one-sigma band and to `1e-4` relative.
- The full objective, data and priors, differs by at most `0.01` in chi-squared: a
  0.1-sigma move along any direction.

Refits of the committed code moved alpha by at most `1.1 %` of its band (He) and
`1.5e-5` relative (Ne), and the objective by at most `1.1e-3` (Ne).

`tests/test_compact.py` re-derives each compact model from its stored parent and
compares the pole structure (count, regions, kept poles and starting point), the
moments, and the worst residuals on the check grid, to `0.01` sigma of the parent's
band. Re-derivations moved CO2's pole strengths by 4.4 % but its residuals by at most
`2.2e-3` sigma.

## Measuring test coverage

The development gate requires **100% statement and branch coverage** across the entire
`optmat` package, including modules a test never imports. `uv run poe check` runs the
default tests and doctests once under coverage; `uv run poe coverage-check` runs this
part of the gate on its own. Even one missing statement or branch fails the gate.

For focused investigations and the longer verification lanes, these reporting commands
measure the same scope without enforcing a threshold:

```sh
uv run poe coverage --randomly-seed=20260930
uv run poe coverage-fits --randomly-seed=20260930
PAPER_LIBRARY=<path> uv run poe coverage-library --randomly-seed=20260930
```

The default command runs the usual tests and doctests. The other two select the
existing `fits` and `library` markers; full refits are deliberately separate because
instrumentation can make them substantially slower. Additional arguments are passed
to pytest, so a test path can narrow a diagnostic run. The seed makes test order
repeatable; use a second seed to check isolation.

The main suite keeps its runtime type checker. Integration tests also exercise
ordinary Python processes, where those development-only checks are absent, to verify
the package's own public validation. Coverage is measured with `pytest-cov`, which
combines the data of the `pytest-xdist` workers; coverage's subprocess patch
(`patch = ["subprocess"]` in `pyproject.toml`) measures those ordinary-Python
processes too, in the same package scope.

Each command writes `.coverage` at the repository root and its reports under
`build/coverage/`: the terminal summary (files with missing lines only), an HTML report
(`html/`, `fits/` or `library/`) and, for the gate, `coverage.xml`, the report CI uploads.
A failed test still makes the command fail. The library command requires
`PAPER_LIBRARY`. Reports and extraction output remain untracked.

The denominator is package code in `src/optmat/`, excluding the coverage tool's existing
default exclusions, including existing `pragma: no cover` annotations. Tests,
extraction scripts, fitting commands and other development tools are outside that
denominator; their verification suites still run. The gate retains the existing
exclusions. Do not add exclusions, remove useful validation or weaken assertions to
meet the threshold. A diagnostic report does not replace the gated default suite.

Statement coverage counts executable statements that ran. Branch coverage counts the
possible exits from decisions that ran: testing a successful validation alone does not
exercise its rejection path. The terminal `Cover` column combines statements and
branch outcomes; JSON reports also give the separate totals. Coverage shows where to
look for missing tests, while assertions establish whether the behaviour is correct.

To combine lanes, run the later ones with `--cov-append`, which adds to the data of the
earlier ones instead of replacing them, then report the total:

```sh
uv run poe coverage
uv run poe coverage-fits --cov-append
PAPER_LIBRARY=<path> uv run poe coverage-library --cov-append
uv run coverage report --fail-under=0
```

Use only lanes run from the same checkout, and identify any omitted lane when reporting
combined coverage. Keep the default-suite figures visible separately.
