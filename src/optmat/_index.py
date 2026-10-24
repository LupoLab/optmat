# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The index of the curated catalogue: which record holds which datasets and formulas.

A search by material, quantity or formula would otherwise read every source record.
The index, ``data/catalogue_index.json``, lists for each record its size, the
materials and quantities of its datasets and the materials of its formulas, so a search
reads only the records that match. It is derived from the records alone (never their
data files) by :func:`build`, and written by ``tools/catalogue_index.py``; a test keeps
the committed file equal to what :func:`build` gives.

The index is checked against the records when first used, at the cost of one ``stat``
per record: a record that is new or whose size differs is read again (its TOML only),
so a search is right even when the index is out of date. If a record's entries then
differ from the index's, or a record has gone, a :class:`CatalogueIndexWarning` says
so. An edit that keeps a record's size (one material id for another of the same
length) is not seen here; the test of the committed index sees it.
"""

import functools
import json
import tomllib
import warnings
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, Final

from optmat.datasets import DATA_DIRECTORY
from optmat.exceptions import CatalogueIndexWarning

INDEX_FILE: Final = DATA_DIRECTORY / "catalogue_index.json"
"""The committed index."""

INDEX_FORMAT: Final = 1
"""The version of the index's layout; an index of another version is not used."""

_SOURCES: Final = DATA_DIRECTORY / "sources"


def entry(path: Path) -> dict[str, Any]:
    """Return one record's index entry, read from its TOML alone.

    Parameters
    ----------
    path : pathlib.Path
        The record, ``data/sources/<source>.toml``.

    Returns
    -------
    dict
        ``size`` (bytes), ``datasets`` (each ``[material, quantity]`` its datasets
        have) and ``formulas`` (each material its formulas are of), in the order
        first met in the record.
    """
    text = path.read_bytes()
    record = tomllib.loads(text.decode("utf-8"))
    pairs = (
        (str(d.get("material", "")), str(d.get("quantity", "")))
        for d in record.get("datasets", [])
    )
    formulas = (str(f.get("material", "")) for f in record.get("formulas", []))
    return {
        "size": len(text),
        "datasets": [list(pair) for pair in dict.fromkeys(pairs)],
        "formulas": list(dict.fromkeys(formulas)),
    }


def build() -> dict[str, Any]:
    """Return the index of the records on disk.

    Returns
    -------
    dict
        ``format`` and ``sources``, each record's :func:`entry` by id in id order.
    """
    paths = sorted(_SOURCES.glob("*.toml"))
    return {"format": INDEX_FORMAT, "sources": {p.stem: entry(p) for p in paths}}


def dumps(index: Mapping[str, Any]) -> str:
    """Return an index as the committed file's text, one record to a line.

    Parameters
    ----------
    index : mapping
        An index from :func:`build`.

    Returns
    -------
    str
        The JSON text, ending in a newline.
    """
    lines = [
        f"{json.dumps(source)}: {json.dumps(value, ensure_ascii=False)}"
        for source, value in index["sources"].items()
    ]
    body = ",\n".join(lines)
    return f'{{\n"format": {index["format"]},\n"sources": {{\n{body}\n}}\n}}\n'


def _committed() -> Mapping[str, Any]:
    """Return the committed index's entries, or none if it is missing or another format.

    Returns
    -------
    mapping
        Each record's entry, by id.
    """
    try:
        index = json.loads(INDEX_FILE.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    if index.get("format") != INDEX_FORMAT:
        return {}
    return index["sources"]


@functools.cache
def index() -> dict[str, dict[str, Any]]:
    """Return each record's entry, the committed index checked against the records.

    Returns
    -------
    dict of str to dict
        The entries by id, in id order.

    Warns
    -----
    CatalogueIndexWarning
        If a record is new, has gone or holds other datasets or formulas than the
        committed index says; the records on disk are used.
    """
    committed = _committed()
    entries: dict[str, dict[str, Any]] = {}
    stale: list[str] = []
    for path in sorted(_SOURCES.glob("*.toml")):
        listed = committed.get(path.stem)
        if listed is None or listed["size"] != path.stat().st_size:
            found = entry(path)
            if listed is None or {**listed, "size": found["size"]} != found:
                stale.append(path.stem)
            listed = found
        entries[path.stem] = dict(listed)
    stale.extend(sorted(set(committed) - set(entries)))
    if stale:
        shown = ", ".join(stale[:5]) + (", ..." if len(stale) > 5 else "")
        msg = (
            f"optmat's catalogue index ({INDEX_FILE.name}) is out of date for "
            f"{len(stale)} source record(s) ({shown}); the records themselves are "
            "used. Regenerate it with `python tools/catalogue_index.py` in a source "
            "checkout."
        )
        warnings.warn(msg, CatalogueIndexWarning, stacklevel=2)
    return entries


def sources_with(
    material: str | None = None, quantity: str | None = None
) -> tuple[str, ...]:
    """Return the records holding a dataset of a material and a quantity.

    Parameters
    ----------
    material : str, optional
        Only datasets of this material id.
    quantity : str, optional
        Only datasets of this quantity.

    Returns
    -------
    tuple of str
        The record ids, in id order.
    """
    return tuple(
        source
        for source, listed in index().items()
        if any(
            (material is None or m == material) and (quantity is None or q == quantity)
            for m, q in listed["datasets"]
        )
    )


def sources_with_formulas(materials: Iterable[str] | None = None) -> tuple[str, ...]:
    """Return the records with formulas, or with formulas of some materials.

    Parameters
    ----------
    materials : iterable of str, optional
        Only records with a formula of one of these material ids.

    Returns
    -------
    tuple of str
        The record ids, in id order.
    """
    wanted = None if materials is None else frozenset(materials)
    return tuple(
        source
        for source, listed in index().items()
        if any(wanted is None or m in wanted for m in listed["formulas"])
    )
