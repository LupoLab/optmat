# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The built catalogue accounts for the records without inventing their meaning."""

import importlib.util
import json
import re
import tomllib
from pathlib import Path

from docutils.core import publish_doctree

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "source_catalogue", ROOT / "docs/_ext/source_catalogue.py"
)
assert SPEC is not None and SPEC.loader is not None
catalogue = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(catalogue)


def test_every_source_and_entry_appears_once(tmp_path):
    records_dir = ROOT / "src/optmat/data/sources"
    records = {
        p.stem: tomllib.loads(p.read_text()) for p in sorted(records_dir.glob("*.toml"))
    }
    assert catalogue.generate_catalogue(records_dir, tmp_path) == tuple(records)
    assert {p.stem for p in tmp_path.glob("*.rst")} == {*records, "index"}
    index = (tmp_path / "index.rst").read_text()
    assert re.findall(r"^   ([a-z0-9_]+)$", index, re.MULTILINE) == list(records)
    for source_id, record in records.items():
        page = (tmp_path / f"{source_id}.rst").read_text()
        # These blocks retain structured metadata exactly, including corrections;
        # code blocks do not reinterpret punctuation as Markdown or mathematics.
        blocks = re.findall(r"\.\. code-block:: text\n\n((?:   .*\n?)+)", page)
        metadata = [
            json.loads("\n".join(s[3:] for s in b.splitlines())) for b in blocks
        ]
        entries = record.get("datasets", []) + record.get("formulas", [])
        assert len(metadata) == len(entries) + bool(record.get("errata"))
        for entry, rendered in zip(entries, metadata, strict=False):
            assert rendered == {
                k: v for k, v in entry.items() if k not in {"notes", "points"}
            }
        for kind in ("datasets", "formulas"):
            count = len(record.get(kind, []))
            label = kind if count != 1 else kind[:-1]
            assert f"{count} {label}" in page
        if record.get("errata"):
            assert metadata[-1] == record["errata"]
        if "doi" not in record:
            assert "No DOI is recorded." in page


def test_generation_is_incremental_and_removes_stale_pages(tmp_path):
    records = tmp_path / "records"
    records.mkdir()
    source = records / "demo_2026.toml"
    source.write_text('citation = "D. Démo (2026)"\nnotes = "Only notes."\n')
    out = tmp_path / "generated"
    catalogue.generate_catalogue(records, out)
    before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in out.iterdir()}
    catalogue.generate_catalogue(records, out)
    assert before == {
        p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in out.iterdir()
    }
    page = (out / "demo_2026.rst").read_text()
    assert "No datasets are recorded." in page
    assert "No formulas are recorded." in page
    assert "Only notes" in page
    source.unlink()
    catalogue.generate_catalogue(records, out)
    assert {p.name for p in out.iterdir()} == {"index.rst"}


def test_record_prose_is_rendered_as_text():
    # Notes must survive Unicode, source notation, and strings resembling directives.
    text = "Démo *not emphasis* | n_2 < 1.\n\n.. warning:: not a directive"
    rendered = catalogue._paragraph(text)
    document = publish_doctree(rendered)
    assert document.astext() == text
    assert not list(
        document.findall(lambda node: node.tagname in {"warning", "emphasis"})
    )
