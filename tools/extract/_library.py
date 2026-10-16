# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Shared helpers of the extraction scripts: read library tables, write data files.

Each ``tools/extract/<source>.py`` defines ``extract() -> dict[str, list[Point]]``
mapping a data file, ``<gas>/<source>__<dataset>.csv`` as the source record's ``file``
names it, to its points, taken from the extracted tables of the curators' paper
library, and a ``main()`` that writes them. The file lives in ``data/values/<family>/``,
the family of the quantity of the dataset that names it (:func:`path`). The
library-marked test re-runs every ``extract()`` and compares with the committed files,
so a file is always what its script makes of the library.

The library is a private collection of the papers, each in a directory named by the
source record's ``library`` key, with its tables extracted under ``tables/`` and the
data files that go with the paper (supporting information, an author's data file, a
database's table) kept byte for byte under ``data/NN/`` (paperextract's ``attach``, with
the file's origin in the paper's ``extraction.json``). It is found at
``$PAPER_LIBRARY``; without it the committed files are the data, and nothing else in
the package needs it.
"""

import csv
import functools
import os
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

from optmat.datasets import data_file

LIBRARY: Final = Path(os.environ.get("PAPER_LIBRARY", "PAPER_LIBRARY-is-not-set"))
"""The paper library (``$PAPER_LIBRARY``); a path that names the variable when unset."""
DATA: Final = Path(__file__).resolve().parents[2] / "src" / "optmat" / "data"
COLUMNS: Final = ("row", "x", "value", "sigma", "printed")


@dataclass(frozen=True)
class Point:
    """One extracted value: row id, the numbers as printed, and the raw cells."""

    row: str
    x: str
    value: str
    sigma: str = ""
    printed: str = ""
    extra: Mapping[str, str] = field(default_factory=dict)
    """The dataset's extra columns (``x_upper``, a per-row condition), as printed."""


def attachment(directory: str, name: str) -> Path:
    """Return a data file attached to a library paper, by its file name."""
    found = sorted((LIBRARY / directory).glob(f"data/*/{name}"))
    if len(found) != 1:
        msg = f"{directory}: {len(found)} attached files named {name!r}, not one"
        raise FileNotFoundError(msg)
    return found[0]


def table(directory: str, table_id: str) -> list[list[str]]:
    """Return a library table's rows (header first) as lists of cell strings."""
    path = LIBRARY / directory / "tables" / f"{table_id}.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


MINUS: Final = "\N{MINUS SIGN}"
_TIMES_POWER = re.compile(r"^\s*([-+]?[\d.]+)\s*\\times\s*10\^\{?([-+]?\d+)\}?\s*$")


def number(cell: str) -> str:
    r"""Return the number printed in a cell, exactly as printed.

    Strips whitespace and a trailing ``\times 10^{k}`` column factor (the caller's
    dataset scale accounts for it), turns a Unicode minus into ``-`` and drops a
    leading ``+``. Raises ``ValueError`` if what remains is not a plain number.
    """
    text = cell.strip().replace(MINUS, "-")
    match = _TIMES_POWER.match(text)
    if match:
        text = match.group(1)
    text = text.removeprefix("+")
    float(text)  # raises ValueError on anything that is not a plain number
    return text


@functools.cache
def _named() -> dict[str, tuple[str, tuple[str, ...]]]:
    """Return the quantity and extra columns of the dataset that names each file."""
    named = {}
    for record in sorted((DATA / "sources").glob("*.toml")):
        with record.open("rb") as handle:
            for entry in tomllib.load(handle).get("datasets", []):
                if "file" in entry:
                    columns = tuple(entry.get("columns", {}))
                    named[entry["file"]] = (entry["quantity"], columns)
    return named


def columns(relative: str) -> tuple[str, ...]:
    """Return a data file's columns: the five common ones, then the declared extras."""
    return COLUMNS + _named()[relative][1]


def path(relative: str) -> Path:
    """Return where a data file lives, ``data/values/<family>/<relative>``.

    The family is that of the quantity of the dataset whose record names the file;
    a file no dataset names is an error (add its dataset to the record first).
    """
    named = _named()
    if relative not in named:
        msg = f"no dataset names the file {relative!r}"
        raise KeyError(msg)
    return data_file(named[relative][0], relative)


def write(files: dict[str, list[Point]]) -> None:
    """Write each data file where :func:`path` puts it."""
    for relative, points in files.items():
        target = path(relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            header = columns(relative)
            writer.writerow(header)
            for point in points:
                if set(point.extra) != set(header[len(COLUMNS) :]):
                    msg = f"{relative}: row {point.row!r} does not fill {header}"
                    raise ValueError(msg)
                writer.writerow(
                    [point.row, point.x, point.value, point.sigma, point.printed]
                    + [point.extra[name] for name in header[len(COLUMNS) :]]
                )
        print(f"wrote {len(points)} points to {target.relative_to(DATA)}")
