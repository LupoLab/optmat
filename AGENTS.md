# AGENTS.md

Working instructions for coding agents (and people) in this repository.

optmat is a library of curated, validated optical properties of materials. Its value is
the trustworthiness of its data: every number traceable to its source, every correction
justified, every model's uncertainty stated. The rules below exist to keep it that way.
The coding standards are conventional; **the data rules and the review register are
not negotiable.**

**Read first:**
- [`CONTRIBUTING.md`](CONTRIBUTING.md): the short guide.
- The manual pages *The curated data* (`docs/data.md`) and *The review register*
  (`docs/register.md`): the long form of the data format and of the register.
- [OEP 0](docs/oeps/oep-0000-process.md), when the work is a major change: the
  enhancement-proposal process (below).

Rules marked ⚙ are enforced by the configuration in `pyproject.toml` and `REUSE.toml`
(ruff, basedpyright, numpydoc, codespell, reuse, pytest); if one of them and the configuration
disagree, the configuration wins and this file is out of date. When a judgement call
concerns physics or the reading of a paper rather than style, ask rather than guess.

---

## Project facts

| | |
|---|---|
| Package | `optmat` (`src/optmat/`) |
| Purpose | Curated optical properties of materials: gases (dispersion, absorption, third-order and Raman response, ionisation), and glasses, crystals and metals (dispersion; fused silica's Kerr and Raman response) |
| Layout | A medium-independent core (`optmat.datasets`, `optmat.vocabulary`, `optmat.materials`, `optmat.conditions`, `optmat.constants`, `optmat.exceptions`) and one subpackage per kind of medium (`optmat.gas`, `optmat.solid`) |
| Data | `src/optmat/data/`: `sources/` (one record per paper), `values/<family>/<material>/` (extracted values), `external/` (third-party files), `materials.toml` (the material registry), `gas/` and `solid/` (model configurations and stored fits) |
| Register | [`REGISTER.md`](REGISTER.md): the review register, cited as `TC-n` |
| Enhancement proposals | `docs/oeps/`: OEPs for major changes; the process is OEP 0 (`oep-0000-process.md`), the form `template.md` |
| Python | 3.12 or later (`requires-python`); development, the lock and the gate on 3.14 |
| Licences | Code, source records and fits MIT; prose and figures CC-BY-4.0; published data (`data/values/`, `data/external/`) remain the publishers', shipped under the third-party data notice (`REUSE.toml`, `LICENSES/LicenseRef-optmat-third-party.txt`) |
| CI | [`.github/workflows/ci.yml`](.github/workflows/ci.yml): the gate below on Linux and macOS, the tests on Windows, and the tests on Python 3.12 with the lowest supported dependencies |
| Type checker | `basedpyright`, pinned exactly in the dev group |

## Commands

| Task | Command |
|---|---|
| Sync environment | `uv sync` |
| Full gate | `uv run poe check` |
| Tests (default lane, on every core) | `uv run poe test` |
| Focused tests, serial | `uv run pytest tests/test_x.py` |
| Tests with the 100% statement and branch coverage gate | `uv run poe coverage-check` |
| Copyright and licence of every file | `uv run poe reuse` |
| Preview provisional register renumbering against main | `uv run poe register-sync` |
| Re-run the extraction against the paper library | `PAPER_LIBRARY=<path> uv run poe test-library` |
| Re-run every stored fit and compare | `uv run poe test-fits` |
| Build the wheel and sdist and install them offline (the release test) | `uv run poe test-release` |
| Fit one gas | `uv run python tools/fit/run.py <gas>` (then `tools/fit/robustness.py <gas>`) |
| Derive a molecule's compact model from its stored fit | `uv run python tools/fit/compact.py <gas>` |
| Fit one gas's hyperpolarizability | `uv run python tools/fit/nonlinear.py <gas>` (then `tools/fit/nonlinear_robustness.py <gas>` and `tools/fit/nonlinear_support_record.py <gas>`) |
| Where the hyperpolarizability fits are supported (the manual's tables) | `uv run python tools/fit/nonlinear_support.py --update-docs` |
| Derive the hyperpolarizability support records (after a refit, or a change to their selection or data) | `uv run python tools/fit/nonlinear_support_record.py [<gas> ...]` |
| Build docs (runs the gallery) | `uv run poe docs` |

Environment variables used by the tools only (the package itself needs none):
- `PAPER_LIBRARY`: the curators' paper library, whose extracted tables the extraction
  scripts read.
- `HITRAN_API_KEY` or `HITRAN_API_KEY_FILE`: a HITRAN API key, for
  `tools/hitran/band_strengths.py`. The key is never printed or committed.
- `RAJ_2018_PROGRAM`: Raj et al.'s published program, for regenerating that source's
  program files.

## Definition of done

`uv run poe check` runs, in order:
1. `ruff check --fix`;
2. `ruff format`;
3. `ruff check`;
4. `basedpyright`;
5. `numpydoc lint`;
6. `codespell`;
7. `reuse lint`: every file states its copyright and licence;
8. `poe coverage-check`: `pytest` (doctests included) on every core
   (`pytest-xdist`, `-n auto --dist loadfile`) under `pytest-cov`, requiring **100%
   statement and branch coverage** across `src/optmat/`, including the workers and
   ordinary-Python subprocesses (coverage's subprocess patch);
9. `sphinx-build -W -n --keep-going` (the gallery executes).

All nine must pass. CI runs exactly this sequence on every push and pull request; if
you change one, change the other in the same commit. There is no other safety net: CI
runs only what `poe check` runs, so the lanes below are still run by hand. Do not
declare work finished while a step fails, and never suppress a check to make it pass.
If a step genuinely cannot pass, say so and why.

The coverage gate uses the existing package scope and exclusions. Tests, extraction
scripts and development tools are outside its denominator. Do not lower the threshold,
add exclusions or remove useful checks to make it pass: tests must establish intended
behaviour with meaningful assertions. `poe coverage`, `coverage-fits` and
`coverage-library` remain diagnostic commands without thresholds; they do not substitute
for the gate. Reports (terminal, XML for CI, HTML) are in `build/coverage/`; the
default tests run once in the gate.

The `poe` test and coverage tasks run in parallel; the library and fits lanes are
distributed test by test (`--dist load`), as each is one file of independent runs.
Tests must therefore not depend on their order or share writable files; write to
`tmp_path`. A plain `pytest` stays serial for focused work.

Also, before merging:
- **A change to extracted data:** `uv run poe test-library` (needs `PAPER_LIBRARY`).
- **A change to packaging, the data files' layout or the public surface:** `uv run poe
  test-release`.
- **A change to a fit, the fitting code or the data a fit uses:**
  - re-run `tools/fit/run.py` for the fits it changes, and for the fits that hold them
    as a phase-matching partner, in the order given on the manual page *Fitting
    dispersion models*;
  - then `tools/fit/robustness.py` for each of them;
  - then `tools/fit/compact.py` for each of them that has a compact model (CO, CO2,
    CH4, N2O, SF6; `data/gas/compact/`), whose freshness test fails until it is
    re-derived;
  - then `uv run poe test-fits`.

  A hyperpolarizability refit (`tools/fit/nonlinear.py`) is followed by
  `tools/fit/nonlinear_robustness.py` and `tools/fit/nonlinear_support_record.py` for
  that gas; the support record's freshness test fails until it is re-derived.
- **Pages whose tables are tested against the stored fits** (the accuracy pages, the
  third-order manual): regenerate their tables and update their text where a finding
  changes.

New functionality is not done until it is **tested**, **documented in the manual** as
well as in docstrings, **recorded in the register** (below) and noted in `CHANGELOG.md`.
Functionality built from an enhancement proposal is also not done until its phase report
is appended to the OEP (below). In your summary, say what you verified against a primary
source and what you could not.

**Proposal-only changes.** A change that only adds or edits OEPs, their index or their
roadmap entries does not need the full gate. Run instead:
- `git diff --check`;
- `uv run codespell` on the changed files, and `uv run poe reuse` if files were added;
- `uv run pytest tests/test_proposals.py tests/test_register.py`;
- `uv run poe docs`.

---

## The review register

[`REGISTER.md`](REGISTER.md) is the public record of how the data were checked, and it
is kept for every change that touches data or models. An entry is required for:

- **every new source record:** what was transcribed (tables, rows, columns, units,
  scale, state, medium, uncertainty and its level), what was checked and against what
  (the table image, the paper's own column, a reproduction of a printed value), and
  anything doubtful;
- **every erratum:** the printed value, the correction, the reason and the evidence
  (the erratum itself also goes in the source record);
- **every fit decision:**
  - a dataset fitted, fixed or used as a validator;
  - a scale freed;
  - a prior chosen;
  - a range set.

  With the numbers that motivated it;
- **every finding:**
  - a disagreement between sources;
  - a formula that does not reproduce its paper;
  - a model that departs from a published one;
- **every open question**, and its resolution when it comes.

How to write one:
- **Numbering:** append at the end with the next free number, `## TC-n <source id or
  subject> (<date>)`. Numbers on `main` are permanent: never renumber, remove or rewrite
  the substance of those entries. Numbers introduced only on a topic branch are
  provisional until integration and may be reassigned together with every citation.
  A later entry corrects or supersedes an established one and cites it; the earlier
  entry's status is updated to point to it.
- **Structure:** bullets in the register's style: `**Where:**`, `**Transcribed:**`,
  `**Checked:**`, `**Doubt:**`, `**Decision:**`, and a closing `**Status:**` line.
  - `open`: needs a check or a decision;
  - `verified (date)`: checked against the primary source;
  - `corrected (date)`: an error was found and corrected;
  - `decided (date)`: a choice was made, with its reason.
- **Citing:** cite entries from code comments, data records and documentation as
  `TC-n`. `tests/test_register.py` checks that every cited entry exists.
- **Content:** facts, not process. The register is read by users deciding whether to
  trust a value; write for them.

**Concurrent work and integration:** numbers written on a topic branch are provisional;
do not reserve guessed ranges or change another branch's entries. Before merging, run
the four-step `poe register-sync` recipe in *Maintaining the review register*
([`docs/developers/register.md`](docs/developers/register.md#concurrent-branches)),
which renumbers the branch's entries after current `main` and resolves an appended-entry
conflict. Never apply a global number replacement after mixing both branches' citations.

## Enhancement proposals (OEPs)

The long form is [OEP 0](docs/oeps/oep-0000-process.md); the index is
[`docs/oeps/index.md`](docs/oeps/index.md).

- **When.** An OEP is needed for:
  - a new medium or subpackage, or a new family of properties or models;
  - a new quantity family, or a change to the data format or the source schema;
  - a changed recommended default or convention;
  - a new fitting methodology;
  - a change to the public surface that propagation codes build on;
  - a new core dependency;
  - process changes.

  It is not needed for new sources, errata, register entries, refits within the existing
  method, new gases or solids added by the recipes below, fixes, documentation, tests,
  refactoring, tooling, or optimisations that leave results unchanged. When unsure, ask.
- **Writing one.**
  - Copy `docs/oeps/template.md` to `oep-NNNN-short-slug.md`, with the next free number.
  - Add it to the index, and to the roadmap's draft section (after the accepted
    proposals).
  - Merge it early as `Draft`, from a `docs/oep-NNNN-…` branch.
  - Remove an informal roadmap item only if the OEP covers it fully. Otherwise link the
    OEP, state the partial coverage and keep the remaining work.
  - Reconcile these entries whenever the scope or status changes.
- **An OEP stands alone.** It may cite:
  - optmat's and solprop's own docs, code, tests and proposals, and commits on their
    main branches;
  - the literature, with DOIs;
  - public reports, standards and data services.

  It never cites the source code of third-party projects, local paths or machines,
  private notes or unmerged branches. Theory is re-derived, in SI and in the conventions
  of the manual. Primary sources are named, with which are in hand. A citation not yet
  curated as a source record is marked and checked before acceptance, never written from
  memory. The OEP compares implementation options and recommends one; it is not a
  line-by-line plan.
- **Authorship and review.** The **Author** is always the person among the maintainers
  who initiated and manages the OEP. AI agents that helped draft it are named, with
  model and version, as **Drafting agent**. Each completed review of Part I is recorded
  under **Reviewers** as `Name (human)` or `Model version (AI)`, and **Review status**
  summarises them. **Every OEP must be reviewed by at least one human and one frontier
  AI model before it is accepted.**
- **Lifecycle.**
  - `Draft`.
  - `Accepted`, by the maintainers: it is linked from the roadmap and Part I is frozen.
  - `In progress`.
  - `Final`: moved to `docs/oeps/archive/` with `git mv`.
  - Or `Deferred`, `Rejected`, `Withdrawn` or `Superseded`.

  Keep the header, the index, the folder and the roadmap consistent;
  `tests/test_proposals.py` checks them.
- **Implementing one.**
  - Work in the OEP's phases, each on its own `feature/oep-NNNN-…` branch, with the full
    gate and the extraction, fit and robustness runs the phase's changes require.
  - The data rules and the register apply as to any other change.
  - Tests, docstrings, the manual page and, where it helps, a gallery example are part
    of the work.
- **Reporting.** Each phase appends its report to Part II of the OEP, in the pull request
  that merges it. The report gives:
  - the delivered paths;
  - the deviations, and why;
  - the decisions, with their register entries where they concern data or models;
  - the register range it wrote;
  - the verification actually run;
  - the sources checked and not checked;
  - the limitations;
  - numbered open questions;
  - a handoff.

  Report failures and gaps plainly.
- **Finishing.** At `Final`:
  - each manual page the OEP produced ends with a link to it;
  - its open questions move to the roadmap, to issues or to `open` register entries;
  - the `CHANGELOG.md` line cites it.

  A final OEP is history: the manual describes current behaviour.

## Data rules

1. **Provenance for every number.** Every dataset has a source record,
   `src/optmat/data/sources/<source>.toml`, that names:
   - the paper (citation and DOI) and where in it the values are (table, page, equation);
   - the technique and the measurement conditions;
   - the state the values refer to and how the source reduced them (ideal or real gas;
     with or without Lorentz–Lorenz);
   - the wavelength medium, and the stated uncertainty with its confidence level;
   - `reproduced_from` when the values are another paper's.
2. **Extracted, not typed.** Tabulated values are written by that source's extraction
   script (`tools/extract/<source>.py`) from the extracted tables of the paper, keeping
   the printed strings exactly. Never hand-edit a data file; re-run its script. Values
   that exist only in the text are entered in the source record as inline `points`,
   with the page.
3. **Look at the page.** Take numbers from the extracted tables, never from a text
   rendering, and check the table image wherever a warning or a structural oddity
   (merged cells, braces, lost headers, swapped captions) touches the numbers.
4. **Corrections are errata, never edits.** A misprint, a swapped caption, a
   misaligned row or a wrong printed coefficient is corrected by an `[[errata]]` entry
   (or, for a formula, a `[[formulas.corrections]]` entry) with the printed value, the
   correction, the reason and the evidence, and a register entry. The loader checks
   that the printed value is what the file holds.
5. **No double counting.** A value reproduced from another paper is marked
   `reproduced_from` and never enters a fit beside its original.
6. **Nothing from memory.** Never transcribe a value into data from memory, a review
   or a secondary source when the primary one is not at hand; record the citation as
   missing instead. A value quoted from memory in code or docs is marked as such and
   gets a register entry.
7. **The vocabulary is the contract.** A new kind of value needs a quantity (with its
   family, units, abscissae, required parameters and columns) in `optmat.vocabulary`,
   documented on the data page, before data of that kind are added. Do not bend an
   existing quantity to hold something it does not mean.
8. **Figures, only when nothing is printed** (TC-699). Prefer, in order: a table or the
   text; a supplement or data repository; the figure. Never digitise what the authors
   publish as numbers.
   - **Vector figures** (the PDF draws the plot) are read by the source's extraction
     script from the drawing: the coordinates of the markers and line vertices, with
     the axes calibrated from the drawing's own tick marks, never by eye.
   - **Raster figures** (an embedded image) are read by the script from the image at
     its native resolution, with the calibration and the rule that picks each curve's
     pixels written in the script, and each value's pixel in its `printed` cell.
   - The dataset sets `digitised = "vector"` or `"raster"`, names the library figure in
     `table`, and states its uncertainty: the source's where stated, combined with the
     digitising uncertainty (the marker or line width against the axis; for an image,
     also the pixel).
   - Check every extraction against every number the source prints for the same
     quantity, and against the extracted values overlaid on the rendered figure; record
     both in the register.
   - Figure values are `use = "check"`. Values read from an image are never fitted (the
     loader rejects it); values read from a vector drawing are `primary` only with the
     maintainer's approval, recorded in the register.

## Adding data or a material

- **A new source:**
  1. The record.
  2. The extraction script.
  3. The data files (written by the script).
  4. Errata with evidence.
  5. A register entry.
  6. `poe test-library` and the gate.

  If the source changes a fit, see *Definition of done*.
- **A new gas:**
  1. Its entry in `optmat.gas.gases`.
  2. Its sources.
  3. A fit configuration (`data/gas/fits/<gas>.toml`) and the stored fit.
  4. Its default model with the reason (`data/gas/defaults.toml`).
  5. A gallery page.
  6. The accuracy page's rows.
  7. Register entries for every choice.
- **A new medium (glasses, crystals):**
  1. A subpackage, `optmat.<medium>`, with its own models and configurations under
     `data/<medium>/`.
  2. Shared needs go into the core only when a second medium actually needs them.
  3. **The core never imports a medium's subpackage** (⚙ `tests/test_package.py`; the
     one exception, a dataset's gas species, is to be generalised first; TC-317).

## Conventions

- **SI at every interface:** m, s, Pa, K, m⁻³; polarizability in C m²/V. Wavelengths
  are vacuum wavelengths; `omega = 2 pi c / lambda_vac` in rad/s. The literature's
  units are converted once, at the loader.
- **Per particle:** a gas's dispersion is stored as its dynamic dipole polarizability,
  and the index at any density follows from Lorentz–Lorenz with the refractivity
  virial term. A measurement is reduced with the density its source meant.
- **Third order:** the second hyperpolarizability is Shelton & Rice's Taylor-series
  gamma, in atomic units.

The manual's conventions pages are the long form.

## Code

- **Formatting and linting** are owned by `ruff` ⚙; do not fight the formatter.
  Suppressions are rule-specific and carry a reason ⚙.
- **Types:** annotate every public signature ⚙; basedpyright must report no new errors.
  Arrays are annotated with jaxtyping on `np.ndarray` and checked by beartype during the
  tests. Use frozen dataclasses for structured data.
- **Modules** declare `__all__`, are organised by domain concept (no `utils`), and keep
  functions small with explicit inputs and outputs.
- **Errors:** raise a specific subclass of `OptmatError`; never return a sentinel where
  raising is correct. Warnings are errors in the tests ⚙: assert expected ones with
  `pytest.warns`. A model warns (`ModelRangeWarning`) outside its stated range.
- **Docstrings:** NumPy style ⚙, validated by numpydoc ⚙, with units and valid ranges
  for every numerical parameter. `Examples` run as doctests ⚙. A physical constant,
  coefficient or citation in code cites its primary source.
- **Comments** explain why: the physics, the source, the non-obvious choice, the
  register entry. Not what the line does.
- **Tests** mirror the package under `tests/`, are deterministic and isolated ⚙
  (random order), use no network, and justify every numerical tolerance in a comment.
  Tests that need the paper library are marked `library` and skipped without it.
- **Documentation:** every feature is described in the manual under `docs/` (Sphinx,
  myst, numpydoc, sphinx-gallery), not only in docstrings. The build runs with
  `-W -n` ⚙.
- **Dependencies:** managed with uv (`uv add`; never `pip install` into the environment).
  Ask before adding one. The floors in `pyproject.toml` are the oldest versions the code
  needs, not the newest available: do not raise one without a reason (an API, a fixed
  bug, a value the tests reproduce), recorded in `CHANGELOG.md`. CI tests the floors.
- **Text files are UTF-8, said explicitly:** every text-mode `open`, `Path.open`,
  `gzip.open`, `read_text` and `write_text` passes `encoding="utf-8"` ⚙
  (`tests/test_encoding.py`); Windows otherwise decodes with its code page.
- **Copyright and licence:** every code file (`.py`, `.jl`) starts with the SPDX header
  ⚙ (`reuse lint`), below any shebang and above the module docstring, as written by
  `uv run reuse annotate --copyright "John C. Travers and contributors" --year 2026
  --license MIT <file>`. Data, prose and figures are covered by `REUSE.toml` (prose and
  figures CC-BY-4.0). Published data are the publishers': values in `data/values/` and
  files in `data/external/` fall under the third-party data notice. A new directory
  under `data/external/` gets its own `REUSE.toml` entry with the third party's
  copyright, and its `PROVENANCE.md` states the publisher's terms, or that it states
  none. **The package is self-contained:** it never downloads data at run time; what a
  model needs is shipped with it.
- **Generated files are never hand-edited;** re-run what writes them:
  - `uv.lock` (uv);
  - the docs build: `docs/_build/`, `docs/auto_gallery/`, `docs/api/generated/`,
    `docs/generated/`, `docs/sg_execution_times.rst`;
  - the curated values, `data/values/` (`tools/extract/<source>.py`), and the
    refractiveindex.info records (`tools/convert_refractiveindex_info.py`);
  - stored fits, robustness results, compact models and support records under
    `data/gas/` and `data/solid/` (`tools/fit/`);
  - test fixtures: `tests/data/public_api.json` (`tools/public_api.py`),
    `golden_models.json` and `golden_raman_identities.json` (`tools/golden_models.py`),
    `solprop_*.json` (`tools/solprop_*.py`), `luna_master.json`
    (`tools/luna_reference.jl`);
  - the manual's tables written by `--update-docs` (`tools/fit/nonlinear_support.py`,
    `tools/two_photon_tables.py`).

## Git and review

- **Branches and worktrees:** all development happens on a topic branch in its own
  worktree. Create them without asking: name the branch for its topic (`data/<source>`,
  `fit/<gas>`, `fix/...`, `feature/...`, `docs/...`; `docs/oep-NNNN-…` for an OEP draft
  and `feature/oep-NNNN-…` for each of its phases) and put the worktree under
  `.claude/worktrees/` (ignored), or where your tool puts its worktrees. Continue in an
  existing branch/worktree when asked to. **The main checkout stays on `main`, clean:**
  other local projects build against it, so never edit files, switch branches or leave
  work in it; its one use is to merge a finished branch. Changes reach GitHub as pull
  requests, which CI checks.
- **Review:** a pull request that changes data or models is reviewed for:
  - its register entries, and its errata with their evidence;
  - its extraction reproducing the files;
  - the fits' changes and their stated reasons.

  A pull request that implements a phase of an OEP is also reviewed against the OEP:
  - its phase report;
  - the deviations it records;
  - the header, the index and the history it updates.
- **Commits:** small, logical, with descriptive messages. Stage deliberately
  (`git add <paths>`).
- **Versions:** never tag, release or change a version by hand; versions come from git
  tags.
- **Coding agents:**
  - never push;
  - never run a destructive or history-rewriting command (`git reset --hard`,
    `git checkout .`, `git restore` over uncommitted changes, `git clean`,
    `git stash drop`, `git branch -D`, `git worktree remove --force`,
    `git commit --amend`, `git rebase`, a force push): these can destroy work that is
    not yours;
  - create branches and worktrees freely (above); `fetch`, `pull` and `merge` into
    `main` only when asked or authorised, and do not ask again for an already
    authorised action.
