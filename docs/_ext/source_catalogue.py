# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Build a source catalogue from recorded metadata, without interpreting the data."""

import json
import re
import tomllib
from collections import defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import quote

from sphinx.application import Sphinx

RECORD_URL = "https://github.com/LupoLab/optmat/blob/main/src/optmat/data/sources/"


def _text(value: str) -> str:
    # Record prose is data, never reStructuredText (or an executable directive).
    return re.sub(r"([\\`*_|<>:.+\-#\[\]()!])", r"\\\1", value)


def _paragraph(value: str) -> str:
    return "\n\n".join(" ".join(p.split()) for p in _text(value).split("\n\n"))


def _literal(value: Any) -> str:
    # Heterogeneous source metadata is intentionally kept as recorded, rather than
    # coerced into a second schema owned by the documentation extension.
    rendered = json.dumps(value, ensure_ascii=False, indent=2)
    return ".. code-block:: text\n\n" + "\n".join(
        "   " + line for line in rendered.splitlines()
    )


def render_record(source_id: str, record: dict[str, Any]) -> str:
    """Render one source's metadata, including records without numerical content."""
    datasets = record.get("datasets", [])
    formulas = record.get("formulas", [])
    materials = sorted({d["material"] for d in datasets + formulas if "material" in d})
    quantities = sorted({d["quantity"] for d in datasets})
    title = source_id
    parts = [
        title + "\n" + "=" * len(title),
        _paragraph(record["citation"]),
        f"`Source record <{RECORD_URL}{quote(source_id)}.toml>`__ · "
        ":doc:`Catalogue </the-data/sources>` · :doc:`Review register </register>`",
    ]
    if record.get("doi"):
        doi = record["doi"]
        parts.append(f"DOI: `{_text(doi)} <https://doi.org/{quote(doi, safe='/')}>`__")
    else:
        parts.append("No DOI is recorded.")
    parts.extend(
        [
            "**Source kind:** " + _text(record.get("kind", "not recorded")),
            "**Recorded technique:** "
            + _paragraph(record.get("technique", "Not recorded.")),
            "**Materials in structured entries:** "
            + (_text(", ".join(materials)) if materials else "None recorded."),
            "**Dataset quantities:** "
            + (_text(", ".join(quantities)) if quantities else "None recorded."),
            f"**Contents:** {len(datasets)} "
            f"dataset{'s' if len(datasets) != 1 else ''}; "
            f"{len(formulas)} formula{'s' if len(formulas) != 1 else ''}.",
            "Curation notes\n--------------",
            _paragraph(record.get("notes", "No curation notes are recorded.")),
        ]
    )
    if record.get("identity"):
        parts.append(
            "**Recorded bibliographic identity status:** "
            + _text(record["identity"])
            + ". This concerns identification of the source; it does not establish "
            "the scientific correctness of its contents."
        )
    parts.extend(_entries(source_id, datasets, formulas))
    parts.append("Corrections\n-----------")
    if record.get("errata"):
        parts.append(_literal(record["errata"]))
    else:
        parts.append("No dataset errata are recorded in this source record.")
    parts.append(
        "Formula corrections, where recorded, appear with their formulas above. "
        "Curation notes may also discuss corrections or unresolved interpretation. "
        "Consult the :doc:`review register </register>` for the checks and decisions."
    )
    return "\n\n".join(parts) + "\n"


def _entries(
    source_id: str, datasets: list[dict[str, Any]], formulas: list[dict[str, Any]]
) -> list[str]:
    parts = []
    parts.append("Datasets\n--------")
    if not datasets:
        parts.append("No datasets are recorded.")
    for dataset in datasets:
        heading = source_id + "/" + dataset["id"]
        parts.append(_text(heading) + "\n" + "~" * len(_text(heading)))
        # Inline points and CSV contents remain in the linked record/data files.
        # All other dataset metadata, including locations, uncertainty, exclusions,
        # and reproduction information, is retained here.
        metadata = {k: v for k, v in dataset.items() if k not in {"points", "notes"}}
        parts.append(_literal(metadata))
        if "points" in dataset:
            parts.append(
                f"{len(dataset['points'])} inline points in the source record."
            )
        if "use" not in dataset:
            parts.append("No use designation is explicitly recorded.")
        if dataset.get("notes"):
            parts.append(_paragraph(dataset["notes"]))
    parts.append("Formulas\n--------")
    if not formulas:
        parts.append("No formulas are recorded.")
    for formula in formulas:
        heading = source_id + "/" + formula["id"]
        parts.extend(
            [
                _text(heading) + "\n" + "~" * len(_text(heading)),
                _literal({k: v for k, v in formula.items() if k != "notes"}),
            ]
        )
        if formula.get("notes"):
            parts.append(_paragraph(formula["notes"]))
    return parts


def generate_catalogue(records_dir: Path, output_dir: Path) -> tuple[str, ...]:
    """Write deterministic pages, changing only files whose content has changed."""
    output_dir.mkdir(parents=True, exist_ok=True)
    records = {
        path.stem: tomllib.loads(path.read_text(encoding="utf-8"))
        for path in sorted(records_dir.glob("*.toml"))
    }
    pages = {
        source_id: render_record(source_id, record)
        for source_id, record in records.items()
    }
    groups: dict[str, list[str]] = defaultdict(list)
    for source_id, record in records.items():
        groups[source_id[0].upper()].append(
            f":doc:`{source_id}` — {_paragraph(record['citation'])}"
        )
    index = [
        "Alphabetical source index\n=========================",
        f"{len(records)} source records, ordered by source identifier. "
        "Each record appears once. See :doc:`/the-data/sources` for how to read "
        "the catalogue and its relation to model selection and the review register.",
        ".. contents:: Jump to a letter\n   :local:\n   :depth: 1",
    ]
    for letter, entries in groups.items():
        index.extend([letter + "\n" + "-" * len(letter), "\n\n".join(entries)])
    index.append(
        ".. toctree::\n   :hidden:\n\n" + "\n".join("   " + s for s in records)
    )
    pages["index"] = "\n\n".join(index) + "\n"
    for name, content in pages.items():
        path = output_dir / (name + ".rst")
        if not path.exists() or path.read_text(encoding="utf-8") != content:
            path.write_text(content, encoding="utf-8")
    # This directory contains only this extension's generated pages.
    for path in output_dir.glob("*.rst"):
        if path.stem not in pages:
            path.unlink()
    return tuple(records)


def _build_catalogue(app: Sphinx) -> None:
    docs = Path(app.srcdir)
    generate_catalogue(
        docs.parent / "src/optmat/data/sources", docs / "generated/sources"
    )


def setup(app: Sphinx) -> dict[str, bool]:
    """Generate before Sphinx discovers documents, on clean and incremental builds."""
    app.connect("builder-inited", _build_catalogue)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
