# The review register

The review register is the record of how optmat's data and models were checked. It is
the part of the package that makes its other claims verifiable: every value in the
curated records can be traced to its page, and the register says what was done with
it once it was there.

## What it records

- **Every transcription:**
  - which tables, rows and columns of a paper were taken;
  - in what unit, at what state and in what wavelength medium;
  - with what stated uncertainty, and at what confidence level;
  - how they were checked: against the table image, against the paper's own
    derived column, by reproducing a printed value from its formula.
- **Every correction** of a printed value, with the evidence. The correction itself
  lives in the source record as an erratum; the register explains it.
- **Every decision** taken in fitting:
  - a dataset fitted, fixed or kept as a validator;
  - a scale left free;
  - a prior chosen;
  - a range set.

  Each comes with the numbers that motivated it.
- **Every finding:**
  - two sources that disagree;
  - a formula that does not reproduce its own paper;
  - a published model that departs from the data.
- **Every question left open**, and its resolution when it comes.

## How to read it

Entries are numbered `TC-n` ("to check"), in the order in which the work was done, and
are permanent once integrated into `main`; branch-only numbers are provisional until
integration. Established entries are never renumbered or removed. The code, the data
records and this manual cite them by number: when a comment, a record's notes or a page
of this manual says `TC-92`, the entry of that number holds the evidence. A later entry
that corrects or supersedes an
earlier one says so. Each entry closes with its status:

| Status | Meaning |
|---|---|
| `open` | Waiting for a check or a decision; the entry says what is open. |
| `verified (date)` | Checked against the primary source. |
| `corrected (date)` | An error was found and corrected. |
| `decided (date)` | A choice was made, for the reason given. Many transcription entries are "accepted with" an integration entry: the records were reviewed as a set and adopted; anything the entry marks open stays open. |

The curation was done by J. C. Travers working with AI agents. Most records were
transcribed and checked by an agent at extraction and reviewed by a second curator, an
agent or J. C. Travers; the decisions on data, fits and conventions are his. "At
extraction" and "a second curator" in the entries refer to this division of the work.

Entries describe the data format of their day. The format was generalised on
2026-09-26 (TC-316), and the package layout settled the same day (TC-317); later entries
use the present names.

## How entries are added

Instructions for authors are in {doc}`/developers/register`. The register itself is
[`REGISTER.md`](https://github.com/LupoLab/optmat/blob/main/REGISTER.md).

(concurrent-branches)=
See [reconciling concurrent branches](developers/register.md#concurrent-branches).

## The register

The entries follow, from TC-1.

```{include} ../REGISTER.md
:start-after: "Status: `open`, `verified (date)`, `corrected (date)`, `decided (date)`."
:heading-offset: 1
```
