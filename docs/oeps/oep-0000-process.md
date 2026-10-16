# OEP 0 — Purpose and process

| | |
|---|---|
| **OEP** | 0 |
| **Title** | Purpose and process |
| **Author** | John Travers |
| **Drafting agent** | Opus 5.5 |
| **Status** | Active |
| **Review status** | Not reviewed |
| **Reviewers** | — |
| **Type** | Process |
| **Created** | 2026-09-29 |
| **Resolution** | Adopted 2026-09-29 |
| **Requires** | — |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | — |
| **Manual pages** | — |

## 1. What an OEP is

An optmat Enhancement Proposal (OEP) is the design document for one major change to
optmat. That might be a new medium, a new family of properties or of models, a new kind
of curated value, a changed default or convention, or a redesigned public interface. The
OEP says:

- what is proposed and why;
- the physics and mathematics, with their references;
- the primary sources the work needs;
- how the change could be implemented;
- how it will be shown to be right.

Once accepted, it is implemented on branches, in phases, and each phase appends a report
to the same document. When the work is complete the OEP is archived, and the manual pages
it produced link to it. Every major feature thus carries the record of how it was
specified, built and verified.

The process follows established practice in the scientific Python community:

- Python's PEPs ([PEP 1](https://peps.python.org/pep-0001/), the process, and
  [PEP 12](https://peps.python.org/pep-0012/), the template);
- NumPy's NEPs ([NEP 0](https://numpy.org/neps/nep-0000.html) and the
  [NEP template](https://numpy.org/neps/nep-template.html)).

It follows the enhancement proposals of the propagation code
[solprop](https://github.com/LupoLab/solprop) most closely. It differs from PEPs and NEPs
in one respect: an OEP is also the record of its own implementation. In a data library,
the decisions made while building and qualifying a feature matter as much as the
proposal. Examples are a source preferred to another, a tolerance, a validator held out
of a fit, and a case refused.

optmat's data are already checked in public, in the review register (`REGISTER.md`,
{doc}`../register`). OEPs do the same for design. They apply to new work from the
adoption of this document. No OEP is written retrospectively for functionality that
already exists.

## 2. When an OEP is needed

An OEP is required for:

- a new medium, or a new subpackage of one (as `optmat.solid` was);
- a new family of properties or of models for a medium (as the Raman responses and the
  ionisation rates of the gases were);
- a new quantity family in `optmat.vocabulary`, or a change to the data format or to the
  schema of the source records ({doc}`../data`);
- a change to a recommended default model, or to a convention of the manual: SI at every
  interface, per-particle storage, Shelton & Rice's gamma ({doc}`../conventions`);
- a new fitting methodology, or a change to how uncertainties are stated;
- a change to the public surface that propagation codes and other users build on, or a
  change that breaks many users' scripts;
- a new dependency of the core;
- a change to how optmat is developed (a *Process* OEP, like this one).

An OEP is not needed for:

- a new source record, its extraction or its errata, or a register entry;
- a refit within the existing methodology;
- a new gas or solid added by the recipes of `AGENTS.md` (*Adding data or a material*);
- bug fixes, documentation, tests, refactoring, tooling and continuous integration;
- dependency updates in the development group;
- small additions to the public interface (a keyword, a helper);
- optimisations that leave results unchanged.

Those follow the working agreement of `AGENTS.md`: an issue or a short agreed plan, then
a pull request, with register entries wherever data or models change.

When it is unclear whether a change needs an OEP, the maintainers decide. They may ask
for one on an issue or a pull request that turns out to be larger than it looked.

## 3. Types

| Type | For |
|---|---|
| Feature | New functionality: a medium, a property, a model family, a quantity, a tool |
| Change | A change to existing behaviour: a default, a convention, an interface, the data format |
| Process | How optmat is developed, reviewed or released |

## 4. Statuses

```text
Draft ──► Accepted ──► In progress ──► Final
  │           │             │
  ├──► Deferred ◄───────────┤
  ├──► Withdrawn            └──► (back to Draft for a substantial revision)
  └──► Rejected
Any status ──► Superseded (by a later OEP)
```

| Status | Meaning | Lives in | On the roadmap |
|---|---|---|---|
| Draft | Written and under discussion; may change freely | `oeps/` | yes, under drafts |
| Accepted | Agreed by the maintainers; waiting to be picked up | `oeps/` | yes |
| In progress | Being implemented; at least one phase has started | `oeps/` | yes |
| Final | Implemented, tested and documented | `oeps/archive/` | no |
| Deferred | Not pursued for now; may return to Draft | `oeps/` | no |
| Rejected | Decided against; the reasons are recorded in the OEP | `oeps/archive/` | no |
| Withdrawn | Abandoned by its author | `oeps/archive/` | no |
| Superseded | Replaced by a later OEP, named in the header | `oeps/archive/` | no |
| Active | A Process OEP in force | `oeps/` | no |

The folder, the [index](index.md) and the roadmap are kept consistent with the status in
the header. `tests/test_proposals.py` checks this.

## 5. Roles

- **The author** is the person among the maintainers who initiated the OEP and manages
  it through drafting and acceptance. They are named in the **Author** field. The author
  is always a person.
- **Drafting agents** are the AI agents, if any, that helped draft or revise the
  proposal. They are named with their model and version in the **Drafting agent** field.
  An agent that joins at a later revision is added with the date of that revision. A
  drafting agent is never the author.
- **Implementers** carry out its phases. They may be the author or not, and may be people
  or coding agents. A phase can be picked up by a different implementer from the one
  before, working from the OEP and its reports alone.
- **Reviewers** review the proposal before acceptance, and each phase's pull request.
  Each review of the proposal is recorded in the **Reviewers** field as `Name (human)` or
  `Model version (AI)`, with entries separated by semicolons. The **Review status** field
  summarises them (§6, *Review*). The author may be one of the human reviewers.
- **The maintainers** accept, reject or defer an OEP, and approve the Final status.

## 6. The lifecycle

### Drafting

1. Copy [the template](template.md) to `docs/oeps/oep-NNNN-short-slug.md`, taking the
   next free number from the [index](index.md). Numbers are never reused, even for a
   withdrawn or rejected OEP.
2. Fill in the header and Part I:
   - the **Author**, and any **Drafting agent**;
   - status `Draft`, review status `Not reviewed` and **Reviewers** `—`.

   Leave Part II empty, and start Part III with the date drafted.
3. Add the OEP to the index table, the index's toctree and the *Draft enhancement
   proposals* section of the {doc}`roadmap <../roadmap>`, after the accepted proposals.
   Review the informal roadmap items in the same change:
   - remove an item only when the OEP covers it in full;
   - for partial coverage, link the OEP, state what it covers and keep the remaining
     work.

   Discussion of a possible future extension does not count as coverage.
4. Open a pull request on a branch named `docs/oep-NNNN-short-slug`. A draft merges
   early, once it is coherent rather than once it is agreed, so that it appears in the
   manual and can be discussed.

Revising a draft is free: edit it in place in further pull requests, and add a line to
Part III for each substantial revision. Whenever its scope changes, update its roadmap
summary and the informal items, and restore work that the revised scope no longer
covers.

### Review

A proposal is reviewed while it is a draft. A review is a reading of the whole of Part I
against the acceptance criteria below. Its findings are resolved in the text or recorded
as open questions.

When a review is complete:
- the reviewer is added to **Reviewers** as `Name (human)` or `Model version (AI)`, for
  example `Opus 5.5 (AI)`;
- a line goes into Part III.

**Review status** is one of:

| Review status | Meaning |
|---|---|
| Not reviewed | No review is recorded |
| In review | At least one review is recorded, but not yet one human and one AI review |
| Reviewed | At least one human and at least one frontier AI model have reviewed it |

**Every OEP must be reviewed by at least one human and one frontier AI model before it
is accepted.** `tests/test_proposals.py` checks two things: that an accepted,
in-progress or final OEP has review status `Reviewed`, and that the status agrees with
the recorded reviewers.

### Acceptance

An OEP is ready to accept when:

- the motivation names a concrete use case;
- the theory is complete: every equation is derived or cited, and its approximations and
  range of validity are stated, in the conventions of {doc}`../conventions`;
- the primary sources are identified, with those in hand and those to be requested, and
  the validators, the independent data held out of any fit, are named;
- the implementation options have been compared and one is recommended, with its effect
  on the vocabulary, the data format and the public surface, and its place in the layout
  (the core or a medium's subpackage);
- the verification plan names analytic limits, reproductions of printed values,
  independent references or invariants, and says what will set each tolerance;
- the documentation plan names the manual page, and a gallery example where appropriate;
- the phases are defined, each independently mergeable;
- the open questions are resolved, or explicitly deferred to a named phase;
- its review status is `Reviewed`: at least one human and one frontier AI model have
  reviewed it (*Review*, above).

The maintainers accept an OEP in a pull request that:
- sets the status to `Accepted`;
- fills in **Resolution** (the date, and the pull request or discussion where it was
  decided);
- adds a line to Part III;
- moves the OEP from *Draft enhancement proposals* to *Accepted enhancement proposals* on
  the roadmap.

A rejection or deferral records the reasons in the OEP, removes its proposal listing from
the roadmap and restores any still-relevant informal work that it had replaced. The same
roadmap update applies on withdrawal. For supersession, reconcile the work with the
replacement OEP's scope.

From acceptance, Part I is frozen. A deviation found during implementation is recorded
in the phase report (§8), not edited into the proposal. A change large enough to
invalidate the proposal is handled in one of two ways:
- the OEP goes back to `Draft` for re-acceptance, and its roadmap listing moves back to
  the draft section;
- or the change is made in a new OEP that supersedes it.

### Implementation

An implementer picks up an accepted OEP in a pull request that:
- sets the status to `In progress`;
- records, in the **Implementation** field, who is implementing it and the branch of the
  first phase.

Each phase is then:

1. implemented on its own branch, `feature/oep-NNNN-short-slug` (with a suffix naming the
   phase when there are several), under the rules of §7;
2. qualified with the full gate, `uv run poe check`, together with the extraction, fit
   and robustness runs that `AGENTS.md` requires for the data and fits it touches;
3. reported: the phase report (§8) is appended to Part II in the same pull request, and
   the **Implementation** field and Part III are updated;
4. reviewed and merged.

Phases may be taken by different implementers, and independent phases may run in
parallel.

### Completion

An OEP becomes `Final` when:

- every phase is merged and reported;
- the manual documents it, and each manual page it produced ends with a link to the OEP;
- its data and models have their register entries;
- the open questions left are moved to the roadmap, to issues or to `open` register
  entries, and named in the last report.

In the pull request that sets `Final`:
- fill in **Manual pages**;
- move the file with `git mv` to `docs/oeps/archive/`;
- move its index row to the archived table, and its toctree entry to the archive index;
- fix the relative links the move breaks;
- remove its roadmap link;
- cite the OEP in the `CHANGELOG.md` line.

A Final OEP is a historical record. It is not updated to follow later changes to the
code; the manual describes current behaviour. Only errata and a "Superseded-by" note are
added afterwards.

## 7. Writing and implementing an OEP

### An OEP stands alone

An OEP is written so that a contributor, or a coding agent, who has never seen the
discussion behind it could implement it from the OEP and optmat's own documentation and
code. It may refer to:

- optmat and solprop: their documentation, modules, tests and enhancement proposals, and
  commits on their main branches;
- the scientific literature, with DOIs where they exist;
- public reports, standards and data services (NIST's databases, HITRAN, CIPM, ISO).

It does not refer to:

- the source code of third-party projects;
- local files, machines or paths;
- private notes, including either project's development archive;
- unmerged branches.

A proposal may well be informed by another implementation. The OEP then re-derives what
it takes, chooses deliberately and cites the literature.

Theory is written in SI units and in the conventions of {doc}`../conventions`, with maths
in MyST dollar notation. A departure from a convention is stated where it occurs. The
implementation section proposes strategies and compares options. It is not a
line-by-line plan, which belongs to the implementer.

A reference that is not yet a curated source record is marked as such in the OEP. It is
checked against the paper before acceptance, never written from memory.

### Implementation rules

- **The data rules hold in full.** Every number comes from a source record. Values are
  written by an extraction script and corrected by errata. Nothing is taken from memory.
  A new kind of value gets its vocabulary quantity before any data of that kind
  (`AGENTS.md`, *Data rules*).
- **The register records the facts.** Every source transcribed, erratum, fit decision,
  finding and open question has its register entry, as for any other change. The phase
  report cites the entries by number, and does not repeat them.
- **The layout.** A medium's work stays in its subpackage. The core takes shared needs
  only when a second medium actually needs them. The core never imports a medium.
- **Tests.** Unit and invariant tests, reproductions of printed values, and comparisons
  with the independent references the OEP names. Every tolerance is justified in a
  comment. The manual's tables are tested against the stored fits where they quote them.
- **Fits.** A phase that changes a fit, or the data a fit uses, re-runs the fits,
  robustness records and fit tests as `AGENTS.md` requires, and says what changed.
- **Downstream.** A change that solprop consumes updates the solprop parity fixtures
  under `tests/data/`, and its report says so.
- **Documentation.** Comments on the physics and the sources, docstrings with units and
  ranges on every public name, and a manual page with the background and a worked
  example. Where it helps, a gallery example.

## 8. The implementation record

Part II of an OEP holds one report per phase, appended when the phase merges. The
reports are how an OEP can be continued by someone else, and how the reasoning behind
the finished code is kept. Each report states:

- **Scope:** the phase, the date, the implementer, the branch and the merge commit.
- **Delivered:** the modules, data, tests and pages added or changed, as paths.
- **Deviations from the proposal**, and why each was necessary.
- **Decisions:** choices made during implementation that the proposal left open. A
  decision about data or a model is also a register entry with status `decided`, and the
  report cites its number.
- **Register entries:** the range written by the phase (TC-a to TC-b).
- **Verification:** the gate run, and on which tree. The extraction, fit and robustness
  runs, and their results. The references, measured errors and reproduced values. What
  was not run.
- **Sources checked:** which constants, coefficients and citations were checked against
  their primary source, and which could not be.
- **Limitations:** cases refused, and ranges or materials not qualified.
- **Open questions**, numbered Q1, Q2, … through the whole OEP, each with its owner or
  destination: the next phase, the roadmap, an issue or a register entry.
- **Handoff:** what the next phase needs to know.

A report is written honestly. A check that failed, a test left weaker than planned, or
an assumption that had to be invented is said plainly.

## 9. Relation to the other records

| Record | Holds |
|---|---|
| OEP | The proposal, and the history of how it was implemented and verified |
| Register (`REGISTER.md`) | The facts about data and models: sources, errata, fit decisions, findings and open questions, cited as TC-n |
| Manual (`docs/`) | How the functionality works and how to use it, as it is now |
| {doc}`Roadmap <../roadmap>` | Open work: accepted and draft OEPs, then ideas and known gaps not fully covered by one |
| `CHANGELOG.md` | What changed in each release, citing the OEP |

The theory of an OEP moves into the manual when it is implemented. The manual page
carries the derivation in its current form, and the OEP remains the record of the
proposal.

## 10. Changing this process

This process is changed by a Process OEP. A small clarification may instead be made by a
pull request that edits this document and adds a line to its history.

## History

- 2026-09-29: Drafted and adopted. Modelled on solprop's SEP 0, with the parts that
  concern a data library: primary sources, the register, fits and validators.
