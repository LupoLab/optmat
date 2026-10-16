# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The counts of curated records on the front pages are today's.

``README.md``, the functionality overview and the changelog state how many source
records, datasets, values, errata and published formulas there are, and how many of
them hold data of solids. When a record is added, update the numbers there.
"""

import re
from pathlib import Path

import pytest

from optmat.datasets import sources
from optmat.materials import material

ROOT = Path(__file__).parents[1]
PAGES = ("README.md", "docs/introduction/overview.md", "CHANGELOG.md")


def _text(page: str) -> str:
    """Return a page with its line breaks and Markdown emphasis removed."""
    text = (ROOT / page).read_text(encoding="utf-8").replace("**", "")
    return re.sub(r"\s+", " ", text)


def _counts() -> dict[str, str]:
    records = sources()
    every = [d for s in records.values() for d in s.datasets]
    solid = [
        d for d in every if material(d.material).kind in {"glass", "crystal", "metal"}
    ]
    return {
        "records": (
            f"{len(records)} source records: {len(every)} datasets and "
            f"{sum(len(d.values) for d in every):,} values, "
            f"{sum(len(s.errata) for s in records.values())} errata and "
            f"{sum(len(s.formulas) for s in records.values())} published formulas"
        ),
        "solids": (
            f"{len({d.source for d in solid})} of the records hold data of glasses, "
            f"crystals and metals: {len(solid)} datasets for "
            f"{len({d.material for d in solid})} materials"
        ),
    }


@pytest.mark.parametrize("page", ["README.md", "docs/introduction/overview.md"])
def test_the_front_pages_counts_are_current(page):
    text = _text(page)
    for statement in _counts().values():
        assert statement in text, (page, statement)


def test_the_changelogs_counts_are_current():
    assert _counts()["records"] in _text("CHANGELOG.md")
