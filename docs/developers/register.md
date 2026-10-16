# Maintaining the review register

Read {doc}`/register` for the meaning of entries and their statuses. This page
describes how to add and integrate entries; it does not change the permanence of
established records.

## How entries are added

Every change to data or models adds its entries: see
[`AGENTS.md`](https://github.com/LupoLab/optmat/blob/main/AGENTS.md) and
[`CONTRIBUTING.md`](https://github.com/LupoLab/optmat/blob/main/CONTRIBUTING.md) in the repository. The register itself is
[`REGISTER.md`](https://github.com/LupoLab/optmat/blob/main/REGISTER.md). The test suite checks that every entry cited
anywhere exists and that every entry has a status.

### Concurrent branches

New entries on a topic branch take the next numbers visible there. Those numbers are
provisional: two branches may independently choose the same ones. `uv run poe
register-sync` maps branch-only numbers to free numbers after the current end of local
`main`. It uses Git's merge base to tell branch-only entries from established ones, so
there is no shared reservation or offset file to keep in sync. Do not reserve guessed
ranges or change another branch's entries.

Before merging a topic branch into `main`:

1. **Commit** the branch's work, including all new files, and preview the mapping:

   ```sh
   uv run poe register-sync --onto main
   ```

2. **Apply** it with `--write`, review the mapping and every changed citation, then
   commit. The helper rewrites tracked UTF-8 text in one pass, including source-record
   metadata and links; check untracked notes and non-UTF-8 references separately.

   ```sh
   uv run poe register-sync --onto main --write
   ```

3. **Merge `main` into the branch.** If its appended register entries conflict,
   preview the resolution, then apply it:

   ```sh
   uv run poe register-sync --resolve
   uv run poe register-sync --resolve --write
   ```

   This keeps the incoming entries and appends the already-renumbered branch entries.
   It refuses overlapping numbers, branch edits to established entries and manually
   edited or unusual conflict files; reconcile those by hand, with care for which
   branch owns each citation.

4. **Resolve** any other conflicts, stage the resolution, and run the full gate on the
   combined tree before completing the merge.

Only one branch integrates into `main` at a time. If `main` advances in the meantime,
repeat the reconciliation and the gate before merging back: an assigned number is not
reserved until it reaches `main`. Once there, the numbers are permanent; later
corrections follow the rule of a later entry. The helper never stages, commits, merges
or fetches, and does not by itself authorise those Git actions. Never apply a global
number replacement after mixing both branches' citations.
