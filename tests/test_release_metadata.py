# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The release records agree: the authors, and the citation's version and date.

``CITATION.cff``, ``.zenodo.json`` and ``pyproject.toml`` each name the authors, and a
release commit moves the changelog's ``Unreleased`` section under its version and date
(``docs/developers/releasing.md``); the citation must then carry the same two. The files
are read as text: the citation's few keys need no YAML parser.
"""

import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RELEASED = re.compile(r"^## (\d+\.\d+\.\d+) \((\d{4}-\d{2}-\d{2})\)\s*$", re.M)
"""A released changelog section: ``## 0.9.1 (2026-10-24)``."""


def _citation() -> str:
    return (ROOT / "CITATION.cff").read_text(encoding="utf-8")


def _scalar(text: str, key: str) -> str | None:
    """Return a top-level scalar of the citation, unquoted, or ``None``."""
    found = re.search(rf"^{key}: *(.+?) *$", text, re.M)
    return None if found is None else found[1].strip("\"'")


def test_every_record_names_the_same_authors():
    """The citation, the archive record and the distribution agree, in order."""
    citation = [
        f"{given} {family}"
        for family, given in re.findall(
            r"- family-names: *(.+?)\n +given-names: *(.+?)\n", _citation()
        )
    ]
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    archive = [
        " ".join(reversed(creator["name"].split(", ")))
        for creator in zenodo["creators"]
    ]
    with (ROOT / "pyproject.toml").open("rb") as handle:
        project = [a["name"] for a in tomllib.load(handle)["project"]["authors"]]
    assert citation
    assert archive == citation
    assert project == citation


def test_the_citation_names_the_latest_release():
    """The citation carries the version and date of the changelog's newest section."""
    released = RELEASED.findall((ROOT / "CHANGELOG.md").read_text(encoding="utf-8"))
    assert released
    citation = _citation()
    assert (
        _scalar(citation, "version"),
        _scalar(citation, "date-released"),
    ) == released[0]
