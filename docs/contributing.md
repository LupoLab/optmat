# Contributing

Corrections, new measurements, new calculations, new materials, better models and
documentation fixes are all welcome.
- [`CONTRIBUTING.md`](https://github.com/LupoLab/optmat/blob/main/CONTRIBUTING.md) is the guide.
- [`AGENTS.md`](https://github.com/LupoLab/optmat/blob/main/AGENTS.md) holds the full rules, for people and coding agents
  alike.
- Issues and pull requests go to the [repository](https://github.com/LupoLab/optmat).

Choose the workflow for your change:

- **Report a problem with a value** by opening an issue: which value, what you expected,
  and what you compared against (a paper with its table or page, an independent
  calculation), at what conditions. A suspected misprint in a paper is just as welcome.
- **Add or correct data** through a source record, an extraction script that writes the
  values exactly as printed, errata with their evidence for any misprint, and an entry
  in the review register ({doc}`register`). The format is on {doc}`data`; the workflow is in {doc}`/developers/data`.
- **Change a model or a fit** by editing its configuration, re-running the fit and its
  leave-one-out robustness record, updating the pages whose tables are tested against
  the stored fits, and recording the decision, with the numbers behind it, in the
  register.

- **Propose a major change** as an {doc}`enhancement proposal <oeps/index>` (OEP):
  - a new medium, family of properties or kind of curated value;
  - a changed default or convention;
  - a redesigned interface.

  An OEP sets out the theory, the primary sources, the options and the verification
  plan. It is reviewed by at least one human and one frontier AI model before it is
  accepted, and it records its implementation phase by phase.

Every change passes the full gate, `uv run poe check`: lint, format, types, docstrings,
spelling, copyright and licence (`reuse lint`), the tests with 100 % statement and branch
coverage of the package, and this manual's build. The coverage scope and additional
verification lanes are described on
{doc}`/developers/testing`. Proposal-only changes use the shorter gate specified
in the working rules.
Prefer topic branches and separate worktrees; agents confirm their creation or switching
with the maintainer. Concurrent register entries are reconciled before integration as
described on {doc}`/developers/register`. Contributions are licensed under the
project's MIT licence.
