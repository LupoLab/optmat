# OEP NNNN — Title

<!--
Copy this file to docs/oeps/oep-NNNN-short-slug.md, with the next free number from the
index, and replace every placeholder. Keep the headings: the parts and their order are
the same in every OEP. Delete the italic guidance as you write each section. The process
is OEP 0 (oep-0000-process.md); read it first.
Add the draft to the OEP index and its toctree, and to the roadmap's Draft enhancement
proposals section after the accepted proposals. Remove informal roadmap items only when
fully covered by this proposal; for partial coverage, link the OEP and retain the
uncovered work explicitly. Keep these entries current when the scope or status changes.
-->

| | |
|---|---|
| **OEP** | NNNN |
| **Title** | A short title, in sentence case |
| **Author** | The person among the maintainers who initiated the OEP and manages it |
| **Drafting agent** | — (AI agents that helped draft it, by model and version; add later ones with the date) |
| **Status** | Draft |
| **Review status** | Not reviewed |
| **Reviewers** | — (each completed review of Part I: `Name (human)` or `Model version (AI)`, separated by `;`) |
| **Type** | Feature, Change or Process |
| **Created** | YYYY-MM-DD |
| **Resolution** | — (at acceptance: the date and the pull request or discussion) |
| **Requires** | — (OEPs this one builds on) |
| **Supersedes** | — |
| **Superseded-by** | — |
| **Implementation** | — (once picked up: the implementer, and each phase's branch and state) |
| **Manual pages** | — (at Final: the manual pages that document the functionality) |

## Part I — Proposal

*Part I is written and revised while the OEP is a Draft, and frozen when it is accepted.
It stands alone. It cites optmat's and solprop's own documentation and code, the
scientific literature and public reports and data services. It never cites third-party
source code, local files or private notes (OEP 0 §7).*

### 1. Abstract

*One paragraph: what is proposed, for what use, and how success will be judged.*

### 2. Motivation and scope

*The concrete use case or problem, and who needs it. What is in scope and, explicitly,
what is not, including tempting adjacent work.*

### 3. Context

*What optmat does today in this area:*
- *links to the manual pages and the modules involved;*
- *the register entries that bear on it, cited by number;*
- *the roadmap item this answers, if any;*
- *the prior OEPs it builds on.*

### 4. Theory

*The physics and mathematics, derived or cited to primary sources:*
- *the equations;*
- *the approximations and their range of validity;*
- *the reference states and the wavelength medium;*
- *the units (SI) and the conventions of conventions.md, with any departure stated
  where it occurs.*

*Every constant or coefficient to be transcribed is cited here. Write maths with dollar
notation, for example $\alpha(\omega) = \sum_k f_k / (\omega_k^2 - \omega^2)$.*

### 5. Primary sources and data

*What the work rests on:*
- *the papers, data services or standards required, and for each whether it is in hand
  or must be requested;*
- *what will be transcribed or extracted, and from which tables;*
- *the vocabulary quantities needed, new or existing;*
- *the validators: independent data held out of any fit, against which the result is
  judged;*
- *the uncertainty expected, and how it will be stated.*

*Mark a citation that is not yet a curated source record.*

### 6. Proposed interface

*How a user would call it:*
- *the public names;*
- *their parameters, with units and ranges;*
- *a short example sketch (fence it as `text` until it runs).*

*How it composes with existing models, datasets and helpers, and what is refused.*

### 7. Implementation strategies

*Options, with their trade-offs, and the recommended one; this is not a line-by-line
plan. Cover as applicable:*

- *where the work belongs: the core or a medium's subpackage;*
- *changes to the vocabulary, the data format or the source schema;*
- *fits: the datasets fitted, fixed or held out, the priors and the ranges;*
- *the expected cost in time and memory;*
- *new dependencies (to be agreed before acceptance).*

### 8. Verification and validation

*How it will be shown to be right:*
- *analytic solutions and limits;*
- *reproductions of printed values;*
- *sum rules and invariants;*
- *independent measurements or calculations held out as validators;*
- *parity with solprop where it consumes the result.*

*For each, the quantity compared and what will set the tolerance.*

### 9. Documentation

*The manual page(s) to add or extend, the docstrings, and whether a gallery example
belongs.*

### 10. Compatibility and impact

*The effect on existing users and results:*
- *changed recommended values or defaults, and the size of the change;*
- *stored fits that must be re-run;*
- *breaking changes to the public interface;*
- *the solprop parity fixtures;*
- *register entries that are superseded.*

### 11. Alternatives considered

*Other designs, models or sources, and why they were not chosen.*

### 12. Phases

*The proposed phases, each independently mergeable with a green gate. For each phase:
its deliverables, its register entries and what it must verify.*

### 13. Open questions

*Questions to settle before acceptance, or deferred to a named phase.*

### 14. References

*Numbered references, with DOIs where they exist. Mark those not yet curated as source
records.*

1. A. Author and B. Author, "Title," *Journal* **volume**, page (year).
   doi:10.xxxx/xxxxx

## Part II — Implementation record

*One report per phase, appended in the pull request that merges the phase (OEP 0 §8).
Copy the skeleton below for each phase.*

### Phase 1 report — short title

- **Scope:** phase, date, implementer, branch, merge commit.
- **Delivered:** modules, data, tests and pages added or changed, as paths.
- **Deviations from the proposal:** each, with why it was necessary.
- **Decisions:** choices the proposal left open, with their register entries where they
  concern data or a model.
- **Register entries:** the range written by this phase.
- **Verification:** the gate and the extraction, fit and robustness runs, and on which
  tree; references, measured errors and reproduced values; what was not run.
- **Sources checked:** constants and citations checked against their primary source,
  and those that could not be.
- **Limitations:** cases refused, ranges or materials not qualified.
- **Open questions:** Q1, Q2, … (numbered through the whole OEP), each with its
  destination.
- **Handoff:** what the next phase needs to know.

## Part III — History

- YYYY-MM-DD: Drafted.
