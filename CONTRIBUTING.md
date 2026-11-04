# Contributing to optmat

Thank you for your interest. Corrections, new measurements, new calculations, new
materials, better models and documentation fixes are all welcome. Because optmat's value
is the trustworthiness of its data, a contribution that touches data or models comes
with its evidence; this guide says what that means.

## Check what already exists

Before opening an issue or a pull request, and before starting work on either, look for
what already covers it or touches it:

- the [issues](https://github.com/LupoLab/optmat/issues) and
  [pull requests](https://github.com/LupoLab/optmat/pulls), open and recently closed;
- the [OEPs](docs/oeps/index.md), since one may already propose the change or be about to
  alter the code;
- for a question about a value, the [review register](REGISTER.md) and the source record:
  the disagreement or the doubt may already be on record;
- for a change that the propagation code [solprop](https://github.com/LupoLab/solprop)
  consumes, its [issues](https://github.com/LupoLab/solprop/issues) and its
  [SEPs](https://solprop.readthedocs.io/en/latest/dev/seps/index.html), since it may be
  waiting for the change or may be broken by it.

```sh
gh issue list --state all --search "line mixing"
gh pr list --state all --search "line mixing"
git grep -il "line mixing" docs/oeps REGISTER.md
gh issue list --repo LupoLab/solprop --state all --search "line mixing"
```

If one already covers your work, add to it rather than opening another. If one only
touches it, such as an open pull request that changes the same files, link it from what
you open, even if you need nothing from it: GitHub shows the reference on both sides, and
the next reader then knows the two are connected. Link across the repositories as
`LupoLab/solprop#285`. For a report of a wrong value, a quick search of the issues is
enough.

## Reporting a problem with a value

Use the [report form](https://github.com/LupoLab/optmat/issues/new?template=report_a_value.yml),
which asks for:
- which value or model (the gas, the wavelength, the function you called, the version:
  `python -c "import optmat; print(optmat.__version__)"`);
- what you expected and what optmat gives;
- what you compared against (a paper with its table or page, an independent
  calculation, an instrument), at what conditions.

A suspected misprint in a paper is just as welcome: say where, and why you think so. If
you have a paper or data that optmat lacks, the
[source form](https://github.com/LupoLab/optmat/issues/new?template=offer_a_source.yml)
is for that.

A wrong value belongs here. A propagation that goes wrong with the right values belongs
with [solprop](https://github.com/LupoLab/solprop/issues); if you cannot tell which it
is, open it in either and it will be moved.

That is all a report needs. The maintainers fill in the rest.

## Writing an issue

Planned work lives in [GitHub issues](https://github.com/LupoLab/optmat/issues). So an
issue for planned work, written by a person or an agent, has to let someone pick up the
task from its text alone. [#6](https://github.com/LupoLab/optmat/issues/6) is an example.

[Check what already exists](#check-what-already-exists) first, and link what you find.
Write one issue per coherent piece of work, not one per finding or per file. A major
change starts as an OEP, which has a tracking issue of its own
([OEP 0](docs/oeps/oep-0000-process.md),
[template](.github/ISSUE_TEMPLATE/oep_tracking.md)).

Use these sections, which the [planned work template](.github/ISSUE_TEMPLATE/planned_work.md)
lays out, leaving out any that do not apply:

- **What happens**, or what is missing, so that someone new to this code can follow it.
- **Why it matters:** what a user or a developer gets wrong or pays for today.
- **Reproduction.** For a bug, the smallest script that shows it, which you have run,
  with its output. For a performance issue, the numbers, the machine, the build and the
  command.
- **Context:** links to what a reader needs, namely the code, the manual, OEPs, source
  records, register entries (`TC-n`), pull requests and other issues. Link code as a
  permalink to a commit (press `y` on a file's page on GitHub), and say which commit you
  checked against.
- **To do:** checkboxes that are actions, not restated symptoms. Name what must not change
  or break. Check that a suggested fix would work before writing it down.
- **Done when:** what someone can check to know it is finished, with a tolerance and a
  machine for a measurement. A criterion should cover every box, and one that no pull
  request can show should say what stands in for it.

Be exact about what you know:

- Say what you ran and what you only read ("read, not run"), and label a suspicion as one.
- Give a measurement its machine, width and date, and a count its scope ("259 errata
  in the 762 source records of 0.9.0").
- Mark a choice that is not the contributor's with **Needs a decision:**, saying who
  decides and what you recommend. A physics judgement, the reading of a paper, a
  recommended value or a changed default, and a changed public signature are the
  maintainers' ([`AGENTS.md`](AGENTS.md)).
- Name a value's source: the source record and the table or page, and the register entry
  that bears on it. A value in an issue follows the same rule as a value in the data:
  nothing from memory.
- Issues are public, so leave out local paths, private notes and nicknames for machines.
  Name hardware by its model. Cite a paper rather than pasting its text or its tables, and
  do not attach its PDF: the publishers' copyright applies.

An issue records work to be done; the register records what has been checked. A finding
about data or a model still gets its register entry when the change is made, and an open
scientific question its `open` entry.

### Labels

Issues are filtered by [label](https://github.com/LupoLab/optmat/labels), so an issue
carries:

- one **type**, or two if it really is both: `bug` (wrong behaviour or a wrong result),
  `data` (a source record, an extraction, an erratum or a curated value), `performance`
  (slower or larger than it should be, with no wrong result), `enhancement` (a new
  feature or a request), `documentation`, `refactor` (no change in behaviour) or `testing`
  (tests, test lanes and CI coverage);
- one **priority**: `priority: high` (wrong values, wrong public statements or a user's
  first contact), `priority: medium` (real work that should not wait long) or
  `priority: low` (worth doing, no hurry);
- an **area** where one fits, such as `area: raman` or `area: solids`;
- `oep` on an OEP's tracking issue, `solprop` if it was found through solprop or solprop
  needs it, and `good first issue` if it needs no physics judgement and no paper library.

If you cannot set labels, suggest them in the text and the maintainers will.

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

1. [Check what already exists](#check-what-already-exists), so that you neither repeat
   nor collide with work in hand. Look especially for open pull requests that change the
   same files and for OEPs that would change the same behaviour.
2. Work on a branch named for its topic (`data/<source>`, `fit/<gas>`, `fix/...`,
   `feature/...`, `docs/...`, `chore/...`). A major change starts as an enhancement
   proposal (OEP):
   see [OEP 0](docs/oeps/oep-0000-process.md) and the
   [template](docs/oeps/template.md). That means a new medium, a new family of
   properties or models, a new kind of curated value, a changed default or convention,
   or a redesigned interface. For anything else larger than a focused fix, open an issue
   first so the approach can be agreed.
3. Read [`AGENTS.md`](AGENTS.md): the data rules, the register, and the coding
   standards (for people and coding agents alike).
4. Keep the change to the task. Add tests with the code, document new functionality in
   the manual as well as in docstrings, add register entries for data and model
   changes, and add a line to `CHANGELOG.md` under `Unreleased`.
5. Run `uv run poe check`; with data or fit changes also `poe test-library` (if you
   have the library) and `poe test-fits`. A change to the documentation alone, including
   a new or edited OEP, runs the shorter checks listed in [`AGENTS.md`](AGENTS.md),
   *Definition of done*.
6. Open a pull request saying what changed, why, and how it was verified, and what you
   assumed, skipped or could not check. Search again first, since work may have appeared
   while you were busy. Link the issue it closes (`Closes #N`) and any issue, pull request
   or OEP it overlaps with or depends on, here or in solprop, and give it a label or two:
   its type and, where one fits, its area, as for an issue. If you cannot set labels, the
   maintainers will.

## Review checklist

A pull request is ready to merge when:

- it links the issues, pull requests and OEPs it closes, overlaps with or depends on;
- the applicable checks pass: the documentation checks for a documentation-only change;
  otherwise `uv run poe check`, and `poe test-fits` where a fit or its data changed;
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
