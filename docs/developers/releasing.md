# Releasing

A release publishes one source distribution and one wheel of `optmat` to PyPI, archives
the version on Zenodo and updates the manual. The version comes from the `v*` git tag.

## Release notes between releases

Pull requests add independent Markdown fragments under `changelog.d/`, named
`<short-slug>.<category>.md`
([format](https://github.com/LupoLab/optmat/blob/main/changelog.d/README.md)), instead of
editing `CHANGELOG.md`. After each push to `main`, the workflow
`.github/workflows/release-notes.yml` creates or refreshes one draft pull request,
**Prepare the next release notes**, on `chore/release-after-<last-version>`. The latest
released version comes from `CITATION.cff`. The draft groups every pending fragment under
`Unreleased` in `CHANGELOG.md` and removes the collected fragments on its branch. Released
sections are copied unchanged.

Edit pending fragments in ordinary pull requests, not in the generated draft. A manual
refresh is available with:

```sh
gh workflow run release-notes.yml --ref main
```

The updater checks that its branch is exactly what it generated before writing another
commit, so manual edits stop it without being overwritten; a closed draft, or one marked
ready, pauses it. It never force-pushes, approves, merges, tags or publishes.

Its pushes come from the workflow's own token, which starts no other workflow, so CI runs
on the draft only once someone pushes a commit to it: the release commit below.

## Making a release

1. Mark the rolling pull request ready, **before editing it**, which freezes the automatic
   refreshes. If more changes reach `main` afterwards, merge `main` into the release
   branch and run `uv run poe changelog-collect`, which appends only the fragments still
   present. With no fragments pending there is no rolling pull request: make a release
   branch from `main` instead.
2. Make the release commit on that branch:
   - rename `Unreleased` in `CHANGELOG.md` to `## X.Y.Z (YYYY-MM-DD)`, above an empty
     `## Unreleased`, and give it a short summary of the release, saying whether any
     recommended value or model changed;
   - if the curated records changed, state their counts there, with the sentence that
     `tests/test_counts.py` keeps current in `README.md`;
   - give `CITATION.cff` the same `version` and `date-released`
     (`tests/test_release_metadata.py` checks that they agree).

   A minor version may break the public surface; every break is listed under
   *Changed (breaking)*.
3. Run `uv run poe test-release`, which builds and installs the distributions offline. CI
   runs on the commit.
4. Merge it once CI passes.
5. Tag the merge on `main` and push the tag on its own:

   ```sh
   git tag -a vX.Y.Z COMMIT -m "optmat X.Y.Z"
   git push origin vX.Y.Z
   ```

   The repository's rules let only administrators create a `v*` tag.
6. Publish a GitHub release from the tag, with the version's changelog section as its
   notes:

   ```sh
   gh release create vX.Y.Z --verify-tag --title "optmat X.Y.Z" --notes-file NOTES.md
   ```

   The release starts `.github/workflows/release.yml`, which builds the distributions,
   checks that the built version matches the tag and that the wheel installs and imports,
   and publishes to PyPI through trusted publishing. Zenodo archives the release with the
   metadata of `.zenodo.json`. The concept DOI in `CITATION.cff` and the README names every
   version, so it does not change.
7. Tell solprop if it is waiting on the release: comment on the issue or SEP that asked
   for the change.

## One-time setup

These are done, and recorded here for a rebuild of the repository's settings:

- On PyPI, a trusted publisher for `LupoLab/optmat`, the workflow `release.yml` and the
  environment `pypi`, and a GitHub environment named `pypi`.
- On Zenodo, the repository switched on under *GitHub*.
- Under GitHub's *Settings → Actions → General → Workflow permissions*, *Allow GitHub
  Actions to create and approve pull requests* enabled, with the default token permissions
  read-only. Only the release-notes job asks for `contents: write` and
  `pull-requests: write`. The workflow runs only from `main`.
