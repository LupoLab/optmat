# Maintaining the manual

The manual uses Sphinx, MyST Markdown, numpydoc and sphinx-gallery. Its five entrance
sections are defined by `docs/index.md`; each section owns a nested toctree. Keep
each authored page in one navigation location and cross-reference it from other
chapters. Enable readers to follow the user guide in order and to enter at a specific
model or property.

## Build and examples

```sh
uv run poe docs
```

The output is `docs/_build/html/`. Warnings and missing references fail the build.
The gallery runs Python scripts from `docs/gallery/` and writes generated pages and
figures under `docs/auto_gallery/`. Edit the scripts and their narratives, not the
generated pages. API pages are generated from the overview in `docs/api/index.md`
and the package docstrings.

Runnable Markdown examples are checked by `tests/test_docs_snippets.py`; installation
and propagation examples are also exercised by the clean-install release tests.
When splitting a page, give its examples their own imports and setup and add it to
the snippet test's page list. Retain old document URLs and fragment identifiers as
overviews or links when moving material.

## Numerical tables

The accuracy and validation pages contain marked tables with freshness assertions
in `tests/test_*page.py` and `tests/solid/test_*page.py`. Keep their start/end markers
when moving a table and update the test and its generator together. Regenerate
tables from stored results; do not edit numbers to match a narrative.

The third-order support tables are regenerated with:

```sh
uv run python tools/fit/nonlinear_support.py --update-docs
```

Removed public comparison tables are retained in
`tests/data/historical_documentation_tables.md`, with their numerical tests.
Ionisation and Raman table tests search the relevant split chapters; the two-photon
and third-order support generators update their marked destination pages directly.

Other table generators and their target pages are identified by their marker names
in the documentation tests and `tools/`. A change in a scientific conclusion belongs
in the register with its evidence. Removing an obsolete public comparison does not
justify removing meaningful regression coverage.

## Source catalogue

`docs/_ext/source_catalogue.py` reads every TOML record in
`src/optmat/data/sources/` before Sphinx discovers documents. It writes an alphabetical
index and one stable page per source identifier under `docs/generated/sources/`.
Those files are ignored by Git, and unchanged pages are not rewritten. Deleted
source records lose their generated page at the next build.

Edit the source record to correct its metadata, following the data and register
rules. Do not edit the generated catalogue. Notes, formula-only and notes-only
records, missing DOIs, Unicode and punctuation are supported without a second
bibliography or new schema. The catalogue reports recorded content; it does not
decide which datasets are scientifically suitable for a fit.

## Reviewing changes

Check the rendered page, equations, tables and figures as well as the source diff.
Test both the sequential reading route and direct links from the introduction,
API and gallery. Run the focused documentation tests and the full {doc}`testing`
gate, and report any unavailable primary-source or numerical verification.
