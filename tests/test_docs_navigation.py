# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The manual has five entrance sections and a single home for each authored page."""

import re
from pathlib import Path

DOCS = Path(__file__).resolve().parents[1] / "docs"
TOCTREE = re.compile(r"^```\{toctree\}\n(.*?)^```", re.MULTILINE | re.DOTALL)


def _entries(page):
    entries = []
    for block in TOCTREE.findall(page.read_text()):
        for raw in block.splitlines():
            line = raw.strip()
            if not line or line.startswith(":"):
                continue
            labelled = re.fullmatch(r"(.+) <([^>]+)>", line)
            label, target = labelled.groups() if labelled else (None, line)
            if "://" not in target:
                entries.append((label, target))
    return entries


def test_the_header_has_the_five_agreed_sections():
    assert _entries(DOCS / "index.md") == [
        ("Introduction", "introduction/index"),
        ("User guide", "user-guide/index"),
        ("The data", "the-data/index"),
        ("Developers guide", "developers/index"),
        ("API", "api/index"),
    ]


def test_every_authored_page_has_one_place_in_the_manual():
    # Generated gallery/API pages have their own trees; check the authored hierarchy
    # without making a test depend on whether a documentation build has run locally.
    pages = {
        p.resolve()
        for p in DOCS.rglob("*.md")
        if p.relative_to(DOCS).parts[0] not in {"_build", "auto_gallery"}
        and "generated" not in p.relative_to(DOCS).parts
    }
    seen = set()

    def visit(page):
        assert page not in seen, f"repeated navigation entry or cycle: {page}"
        seen.add(page)
        for _, target in _entries(page):
            child = (page.parent / f"{target}.md").resolve()
            if any(
                child.is_relative_to(DOCS / root)
                for root in ("auto_gallery", "generated")
            ):
                continue
            assert child in pages, f"missing authored page: {page} -> {target}"
            visit(child)

    visit(DOCS / "index.md")
    assert pages == seen, f"pages outside the navigation: {pages - seen}"
