# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Read a sheet of an Excel workbook (``.xlsx``) with the standard library only.

An ``.xlsx`` file is a zip of XML parts: ``xl/workbook.xml`` names the sheets and
``xl/_rels/workbook.xml.rels`` maps each to its ``xl/worksheets/sheetN.xml``; text cells
(``t="s"``) index ``xl/sharedStrings.xml``. A numeric cell stores its value as Excel's
serialisation of a binary double (``8.1696000000000009``); :func:`sheet` writes it back
in the shortest decimal form that reads as the same double (``8.1696``), so nothing is
rounded. How the workbook *displays* a number (a cell's number format, such as ``0.00``)
is presentation and is not applied. Formulas are refused: a cell must hold a value.
"""

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Final

_MAIN: Final = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_RELATIONSHIP: Final = (
    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
)
_COLUMN: Final = re.compile(r"^([A-Z]+)\d+$")


def _column(reference: str) -> int:
    """Return the zero-based column of a cell reference such as ``AB12``."""
    match = _COLUMN.match(reference)
    if match is None:
        msg = f"not a cell reference: {reference!r}"
        raise ValueError(msg)
    index = 0
    for letter in match.group(1):
        index = 26 * index + ord(letter) - ord("A") + 1
    return index - 1


def shortest(text: str) -> str:
    """Return the shortest decimal string that reads as the same double as ``text``.

    Python's ``repr`` of a float is that string (``30.0``, ``0.994``); a whole number
    loses its ``.0`` (``30``), as the workbook shows it.
    """
    value = float(text)
    return repr(value).removesuffix(".0")


def sheet(path: Path, name: str) -> list[list[str]]:
    """Return the rows of the sheet ``name``, each a list of cell strings.

    Rows are in sheet order, including empty ones (as ``[]``); missing cells within a
    row are ``""``. Text cells give their string, numeric cells :func:`shortest` of
    the stored value.
    """
    with zipfile.ZipFile(path) as book:
        strings = [
            "".join(t.text or "" for t in item.iter(f"{_MAIN}t"))
            for item in ET.fromstring(book.read("xl/sharedStrings.xml"))
        ]
        targets = {
            rel.get("Id"): rel.get("Target")
            for rel in ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(book.read("xl/workbook.xml"))
        found = [s for s in workbook.iter(f"{_MAIN}sheet") if s.get("name") == name]
        if len(found) != 1:
            msg = f"{path.name}: no single sheet named {name!r}"
            raise KeyError(msg)
        target = targets[found[0].get(_RELATIONSHIP)]
        root = ET.fromstring(book.read(f"xl/{target}"))
    rows: list[list[str]] = []
    for row in root.iter(f"{_MAIN}row"):
        cells: dict[int, str] = {}
        for cell in row.iter(f"{_MAIN}c"):
            if cell.find(f"{_MAIN}f") is not None:
                msg = f"{name}!{cell.get('r')}: a formula, not a value"
                raise ValueError(msg)
            value = cell.find(f"{_MAIN}v")
            if value is None or value.text is None:
                continue
            kind = cell.get("t", "n")
            column = _column(cell.get("r", ""))
            if kind == "s":
                cells[column] = strings[int(value.text)]
            elif kind == "n":
                cells[column] = shortest(value.text)
            else:
                msg = f"{name}!{cell.get('r')}: unexpected cell type {kind!r}"
                raise ValueError(msg)
        rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
    return rows
