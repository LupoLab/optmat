# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Keep the enhancement proposals consistent with their index, folders and roadmap.

Each OEP under ``docs/oeps/`` carries its status in its header table. The status
decides where the OEP lives (the open folder or ``archive/``), how the index lists it
and whether the roadmap links it (OEP 0, ``docs/oeps/oep-0000-process.md``, section 4).

These tests check:

- that the four agree;
- that every header is complete;
- that the review status agrees with the recorded reviewers;
- that nothing is accepted without a human and an AI review (OEP 0 section 6,
  *Review*);
- that no OEP refers to local files or to the private development archive, because an
  OEP stands alone (OEP 0 section 7).

What an OEP says is for review; only its structure is checked here.
"""

import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "docs"
OEPS = DOCS / "oeps"
ARCHIVE = OEPS / "archive"
ROADMAP = DOCS / "roadmap.md"

FIELDS = (
    "OEP",
    "Title",
    "Author",
    "Drafting agent",
    "Status",
    "Review status",
    "Reviewers",
    "Type",
    "Created",
    "Resolution",
    "Requires",
    "Supersedes",
    "Superseded-by",
    "Implementation",
    "Manual pages",
)
OPEN_STATUSES = frozenset({"Draft", "Accepted", "In progress", "Deferred", "Active"})
ARCHIVED_STATUSES = frozenset({"Final", "Rejected", "Withdrawn", "Superseded"})
TYPES = frozenset({"Feature", "Change", "Process"})
REVIEW_STATUSES = frozenset({"Not reviewed", "In review", "Reviewed"})
# Statuses a proposal can only reach after review: acceptance and everything after it.
REVIEWED_STATUSES = frozenset({"Accepted", "In progress", "Final"})
PARTS = (
    "## Part I — Proposal",
    "## Part II — Implementation record",
    "## Part III — History",
)

FILENAME = re.compile(r"oep-(\d{4})-[a-z0-9]+(?:-[a-z0-9]+)*\.md")
TITLE = re.compile(r"\A# OEP (\S+) — (.+)$", re.MULTILINE)
FIELD = re.compile(r"^\| \*\*(.+?)\*\* \| (.*?) \|$", re.MULTILINE)
DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
REVIEWER = re.compile(r"(.+?) \((human|AI)\)")
INDEX_ROW = re.compile(
    r"^\| \[(\d+)\]\(((?:archive/)?oep-[^)]+\.md)\) \| (.+?) \| (.+?) \| (.+?) \|$",
    re.MULTILINE,
)
ROADMAP_LINK = re.compile(r"\(oeps/((?:archive/)?oep-\d{4}-[^)#]+\.md)")
# Absolute or home-relative paths, and paths into the private development archive,
# which lives in dev/ beside the public tree and is never published.
LOCAL_REFERENCE = re.compile(r"/Users/|/home/|file://|~/|[A-Za-z]:\\|(?<![\w/.-])dev/")


def header(page: Path) -> dict[str, str]:
    """Return the fields of an OEP's header table."""
    return dict(FIELD.findall(page.read_text(encoding="utf-8")))


def reviewers(fields: dict[str, str]) -> list[tuple[str, str]]:
    """Return a header's reviewers as ``(name, kind)`` pairs; kind is human or AI."""
    if fields["Reviewers"] == "—":
        return []
    entries = [entry.strip() for entry in fields["Reviewers"].split(";")]
    matches = [REVIEWER.fullmatch(entry) for entry in entries]
    assert all(matches), f"write each as 'Name (human)' or 'Model (AI)': {entries}"
    return [(match[1], match[2]) for match in matches if match]


def relative(page: Path) -> str:
    """Return an OEP's path relative to ``docs/oeps``, as the index links it."""
    return page.relative_to(OEPS).as_posix()


OEP_FILES = sorted([*OEPS.glob("oep-*.md"), *ARCHIVE.glob("oep-*.md")])
IDS = [relative(page) for page in OEP_FILES]


def test_there_is_a_process_oep():
    assert OEPS / "oep-0000-process.md" in OEP_FILES


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_oep_filename_title_and_header_agree(page: Path):
    match = FILENAME.fullmatch(page.name)
    assert match, f"{page.name} is not named oep-NNNN-short-slug.md"
    title = TITLE.search(page.read_text(encoding="utf-8"))
    assert title, f"{page.name} does not start with '# OEP N — Title'"
    fields = header(page)
    assert int(title[1]) == int(match[1]) == int(fields["OEP"])
    assert fields["Title"] == title[2]


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_oep_header_is_complete(page: Path):
    fields = header(page)
    assert [field for field in FIELDS if field not in fields] == []
    assert fields["Status"] in OPEN_STATUSES | ARCHIVED_STATUSES
    assert fields["Type"] in TYPES
    assert DATE.fullmatch(fields["Created"])
    if fields["Status"] == "Superseded":
        assert fields["Superseded-by"] != "—"
    if fields["Status"] == "Active":
        assert fields["Type"] == "Process"
    assert fields["Author"] not in ("", "—"), "every OEP has a human author"
    assert fields["Drafting agent"] != ""
    assert fields["Review status"] in REVIEW_STATUSES


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_the_review_status_agrees_with_the_reviewers(page: Path):
    fields = header(page)
    kinds = {kind for _, kind in reviewers(fields)}
    expected = (
        "Not reviewed"
        if not kinds
        else "Reviewed"
        if kinds == {"human", "AI"}
        else "In review"
    )
    assert fields["Review status"] == expected


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_nothing_is_accepted_without_a_human_and_an_ai_review(page: Path):
    fields = header(page)
    if fields["Status"] in REVIEWED_STATUSES:
        assert fields["Review status"] == "Reviewed", (
            f"a {fields['Status']} OEP needs a human and a frontier AI review"
        )


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_oep_lives_in_the_folder_its_status_names(page: Path):
    status = header(page)["Status"]
    expected = ARCHIVE if status in ARCHIVED_STATUSES else OEPS
    assert page.parent == expected, f"a {status} OEP belongs in {expected.name}/"


@pytest.mark.parametrize("page", OEP_FILES, ids=IDS)
def test_a_proposal_has_its_three_parts_in_order(page: Path):
    if header(page)["Type"] == "Process":
        pytest.skip("a process OEP is structured by its subject")
    text = page.read_text(encoding="utf-8")
    positions = [text.find(part) for part in PARTS]
    assert -1 not in positions, f"missing one of {PARTS}"
    assert positions == sorted(positions)


def test_oep_numbers_are_unique():
    numbers = [int(header(page)["OEP"]) for page in OEP_FILES]
    assert len(numbers) == len(set(numbers))


def test_the_index_lists_every_oep_with_its_header():
    rows = {
        path: (int(number), title, status, kind)
        for number, path, title, status, kind in INDEX_ROW.findall(
            (OEPS / "index.md").read_text(encoding="utf-8")
        )
    }
    assert sorted(rows) == sorted(IDS)
    for page in OEP_FILES:
        fields = header(page)
        assert rows[relative(page)] == (
            int(fields["OEP"]),
            fields["Title"],
            fields["Status"],
            fields["Type"],
        ), relative(page)


def test_the_index_toctree_lists_every_open_oep():
    text = (OEPS / "index.md").read_text(encoding="utf-8")
    toctree = text[text.index("```{toctree}") :]
    entries = {line.strip() for line in toctree.splitlines()}
    expected = {page.stem for page in OEP_FILES if page.parent == OEPS}
    assert expected <= entries
    assert {"template", "archive/index"} <= entries


@pytest.mark.parametrize(
    ("section", "statuses"),
    [
        ("Accepted enhancement proposals", {"Accepted", "In progress"}),
        ("Draft enhancement proposals", {"Draft"}),
    ],
)
def test_the_roadmap_links_exactly_the_oeps_with_each_status(
    section: str, statuses: set[str]
):
    text = ROADMAP.read_text(encoding="utf-8")
    accepted = text.index("\n## Accepted enhancement proposals\n")
    drafts = text.index("\n## Draft enhancement proposals\n")
    assert accepted < drafts
    start = text.index(f"\n## {section}\n")
    end = text.find("\n## ", start + 1)
    linked = set(ROADMAP_LINK.findall(text[start:end]))
    expected = {
        relative(page) for page in OEP_FILES if header(page)["Status"] in statuses
    }
    assert linked == expected


@pytest.mark.parametrize(
    "page",
    [*OEP_FILES, OEPS / "index.md", OEPS / "template.md", ARCHIVE / "index.md"],
    ids=[*IDS, "index.md", "template.md", "archive/index.md"],
)
def test_an_oep_page_refers_to_no_local_files(page: Path):
    assert LOCAL_REFERENCE.findall(page.read_text(encoding="utf-8")) == []


def test_the_template_carries_every_header_field_and_part():
    template = OEPS / "template.md"
    text = template.read_text(encoding="utf-8")
    assert [field for field in FIELDS if field not in header(template)] == []
    assert [part for part in PARTS if part not in text] == []
    assert header(template)["Status"] == "Draft"
    assert header(template)["Review status"] == "Not reviewed"
