# Contributing to optmat

Thank you for your interest. Corrections, new measurements, new calculations, new
materials, better models and documentation fixes are all welcome. Because optmat's value
is the trustworthiness of its data, a contribution that touches data or models comes
with its evidence; this guide says what that means.

## Reporting a problem with a value

Open an issue saying:
- which value or model (the gas, the wavelength, the function you called, the version:
  `python -c "import optmat; print(optmat.__version__)"`);
- what you expected and what optmat gives;
- what you compared against (a paper with its table or page, an independent
  calculation, an instrument), at what conditions.

A suspected misprint in a paper is just as welcome: say where, and why you think so.

## Setting up

optmat supports Python 3.12 and later; development uses [uv](https://docs.astral.sh/uv/)
and Python 3.14:

```sh
git clone https://github.com/LupoLab/optmat.git
cd optmat
uv sync --locked      # the package, editable, with the development tools
uv run poe check      # lint, format, types, docstrings, spelling, licences, tests, docs
```

## Measuring test coverage

The full gate requires 100% statement and branch coverage of `src/optmat/`.
See [Development checks](docs/developers/testing.md) for the gated and diagnostic
commands, coverage scope, reports and additional verification lanes.

## Branches and concurrent register entries

Work on a topic branch, ideally in its own worktree, and open a pull request; CI runs the
full gate on it. Agents create their branch and worktree without asking and never work
in the main checkout.

Register numbers become permanent when they reach `main`; numbers written on a branch
are provisional. Before merging, renumber the branch's entries with `uv run poe
register-sync`, following the four steps in
[Maintaining the review register](docs/developers/register.md#concurrent-branches).

## Adding or correcting data

The [data-format reference](docs/data.md) and [quantity vocabulary](docs/the-data/quantities.md)
describe the records. [Adding sources and materials](docs/developers/data.md) gives
the workflow:

1. **The source record.** One file per paper, `src/optmat/data/sources/<source>.toml`,
   with the citation and DOI, the technique, the conditions, the reference state and
   how the paper reduced to it, the stated uncertainties, and one entry per dataset.
2. **The values, as printed.** Tabulated values are written by an extraction script,
   `tools/extract/<source>.py`, from the table as the paper prints it; values found
   only in the text go in the record with their page. Printed values are never edited.
3. **Errata.** A misprint is corrected by an erratum in the record (the printed value,
   the correction, the reason and the evidence); the loader checks the printed value.
4. **The register.** Add an entry to [`REGISTER.md`](REGISTER.md): what you transcribed,
   what you checked and against what, what is doubtful. Cite it as `TC-n` where the
   change is made.

The extraction scripts read a private library of the papers (`$PAPER_LIBRARY`). If you
do not have the paper in that form, say so in the pull request: one of the maintainers
re-runs the extraction from the published table.

A value that depends on density or temperature must say which state it refers to and
how its source reduced it. Getting that wrong is behind several of the errors we have
found in published formulas.

## Changing a model or a fit

Fits are configured in `src/optmat/data/gas/fits/` (and `nonlinear/` for the third-order
response) and stored under `data/gas/models/`. A change to a fit is re-run with the
tools in `tools/fit/`, its leave-one-out robustness record regenerated, and the pages
whose tables are checked against the stored fits updated; see
[`AGENTS.md`](AGENTS.md), *Definition of done*. Record the decision and the numbers that
motivated it in the register.

## Making a change

1. Work on a branch named for its topic (`data/<source>`, `fit/<gas>`, `fix/...`,
   `feature/...`, `docs/...`). A major change starts as an enhancement proposal (OEP):
   see [OEP 0](docs/oeps/oep-0000-process.md) and the
   [template](docs/oeps/template.md). That means a new medium, a new family of
   properties or models, a new kind of curated value, a changed default or convention,
   or a redesigned interface. For anything else larger than a focused fix, open an issue
   first so the approach can be agreed.
2. Read [`AGENTS.md`](AGENTS.md): the data rules, the register, and the coding
   standards (for people and coding agents alike).
3. Keep the change to the task. Add tests with the code, document new functionality in
   the manual as well as in docstrings, add register entries for data and model
   changes, and add a line to `CHANGELOG.md` under `Unreleased`.
4. Run `uv run poe check`; with data or fit changes also `poe test-library` (if you
   have the library) and `poe test-fits`. A change that only adds or edits OEPs runs
   the shorter checks listed in [`AGENTS.md`](AGENTS.md), *Definition of done*.
5. Open a pull request saying what changed, why, and how it was verified, and what you
   assumed, skipped or could not check.

## Review checklist

A pull request is ready to merge when:

- `uv run poe check` passes, and `poe test-fits` where a fit or its data changed;
- every new or changed value has its source record, its register entry, and, for a
  correction, its erratum with evidence;
- no value was typed from memory or a secondary source where the primary one was
  available;
- a changed fit states why, with the numbers, and its robustness record and tested
  documentation tables are regenerated;
- the manual, docstrings and changelog line are present;
- a pull request that implements a phase of an OEP appends the phase report to the OEP
  and updates its header, index entry and history. The one that completes it archives
  the OEP and links it from the manual pages it produced;
- an OEP is accepted only after at least one human and one frontier AI model have
  reviewed it.

## Licence

By contributing you agree that your contributions are licensed under the project's
[MIT licence](LICENSE). Published data you add (values transcribed from a source, files
from a data service) remain their publishers' and are included under the
[third-party data notice](LICENSES/LicenseRef-optmat-third-party.txt): record the source,
and its terms where it states them, as the data rules require.
