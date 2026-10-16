# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Jhanwar & Meath (1982): the CO and CO2 constrained DOSDs and dipole sums.

The tables print powers of ten in parentheses (``13.79(3)`` = 13.79e3, ``97.27(-1)``
= 9.727) and a reference or a footnote mark after a value (``0.414[17]``,
``0.001550 a)``, ``13.09b,c)``); numbers are written in e-notation with the digits
unchanged, the cell kept in ``printed``. Images checked against the cells
(2026-09-27); the scan is poor in places (see below).

- Table 1 (``tbl_025ed8b378b2``, cells file; its caption stands in the text above it
  and is not attached, structural): the CO DOSD integrated over 24 energy regions
  (eV), "Recommended values" only (the initial-data column holds several references
  per cell). Row = the region, x and x_upper its edges (``10^3`` -> 1000, ``10^4`` ->
  10000, ``∞`` -> inf).
- Table 2 (``tbl_abfee4bd8e42``, cells file): the CO2 DOSD, "UV and higher energy"
  and "total (including IR)" columns. The "infrared" row (total only, the four IR
  bands' summed strength) goes to its own file, label abscissa. The last region's
  lower edge is extracted as ``10^2`` but the scan's superscript is illegible and the
  regions are contiguous (``10^3 - 10^4`` precedes it; the CO table ends ``10^4 -
  ∞``): it is read as 10^4 (structural).
- Table 3 (``tbl_61308cacb9a8``, CSV) and Table 4 (``tbl_b9eb3fd4a103``, CSV): the
  dipole sums of CO and CO2 (UV and higher energy; for CO2 also the total, printed as
  a second line in each cell, which the extraction joins to the first), and the
  literature values, each with a footnote mark naming its source. Table 4's row
  printed "S_{-5}" in the extraction is S_{-8} (poor scan: the caption and text list
  k = -6, -8, -10, and the pseudo-DOSD reproduces the value as S_-8). In the S_-2 rows
  a pair "X and Y" is the Kirkwood-Unsold result of Mulder et al. (footnote e,
  misprinted "c)" in both tables' cells or footnote lists); the text (p. 7) says the
  Kirkwood value for CO agrees very well with 13.08, so the second of each pair is
  the Kirkwood and the first the Unsold value.
- Table 5 (``tbl_9653e9a5f5a6``): the ten-term pseudo-DOSDs are formulas in the
  record (this script asserts that the record's terms equal the table's cells).
"""

import re
import tomllib
from pathlib import Path

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Jhanwar_1982_DipoleOscillatorStrength"
TABLE_1 = "tbl_025ed8b378b2"
TABLE_2 = "tbl_abfee4bd8e42"
TABLE_3 = "tbl_61308cacb9a8"
TABLE_4 = "tbl_b9eb3fd4a103"
TABLE_5 = "tbl_9653e9a5f5a6"
RECORD = Path(__file__).resolve().parents[2] / (
    "src/optmat/data/sources/jhanwar_1982.toml"
)
ORDERS = ["2", "1", "-1", "-2", "-3", "-4", "-6", "-8", "-10"]
# Literature sources by footnote mark, per table.
CO_MARKS = {
    "a": "mulder_1981",
    "b": "parker_1976",
    "c": "watson_1936",
    "d": "liebmann_1971",
}
CO2_MARKS = {
    "a": "mulder_1981",
    "b": "pack_1974",
    "c": "watson_1936",
    "d": "morrison_1977",
}
_NUMBER = re.compile(r"(?P<mantissa>\d+\.\d+)\s*(?:\((?P<exponent>-?\d+)\))?")
_RECOMMENDED = re.compile(r"\d+\.\d{2}(?:\(-?\d+\))?")
_TOKEN = re.compile(
    r"(?P<first>\d+\.\d+(?:\s*\(-?\d+\))?)"
    r"(?:\s*and\s*(?P<second>\d+\.\d+(?:\s*\(-?\d+\))?))?"
    r"\s*(?P<marks>[a-e](?:,[a-e])*)\)"
)


def _value(text: str) -> str:
    """Return a ``mantissa(exponent)`` number in e-notation, digits unchanged."""
    found = _NUMBER.fullmatch(text.strip())
    if found is None:
        msg = f"unreadable number {text!r}"
        raise ValueError(msg)
    mantissa, exponent = found["mantissa"], found["exponent"]
    return number(f"{mantissa}e{exponent}" if exponent else mantissa)


# Printed powers of ten and infinity as region edges; "10^2" is Table 2's illegible
# last lower edge, 10^4 (module docstring).
_EDGES = {"10^3": "1000", "10^4": "10000", "10^2": "10000", "∞": "inf"}


def _bound(text: str) -> str:
    text = text.strip().rstrip("-").strip()
    return _EDGES[text] if text in _EDGES else number(text)


def _dosd(table_id: str, columns: dict[int, str]) -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, table_id)[1:]:
        if cells[0] == "infrared":
            value = cells[3].split("a)")[0]
            point = Point("infrared", "infrared", _value(value), "", " | ".join(cells))
            files["CO2/jhanwar_1982__infrared_CO2.csv"] = [point]
            continue
        if not cells[1]:
            continue  # the sub-header row of Table 2
        # Table 2's last lower edge: "10^2" in the extraction, 10^4 (docstring).
        low, high = _bound(cells[0]), _bound(cells[1])
        row = f"{low}-{high}"
        for column, path in columns.items():
            text = re.sub(r"\s*(\[\d+\]|\(M\))$", "", cells[column])
            point = Point(
                row, low, _value(text), "", " | ".join(cells), {"x_upper": high}
            )
            files.setdefault(path, []).append(point)
    return files


def _tokens(cell: str) -> list[tuple[str, list[str], str]]:
    """Return (marks, values, text) per literature entry of a cell."""
    return [
        (
            found["marks"],
            [_value(found["first"])]
            + ([_value(found["second"])] if found["second"] else []),
            found[0],
        )
        for found in _TOKEN.finditer(cell)
    ]


def _co_sums() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = table(DIRECTORY, TABLE_3)[1:]
    for order, cells in zip(ORDERS, rows, strict=True):
        printed = " | ".join(cells)
        point = Point(order, order, _value(cells[1]), "", printed)
        files.setdefault("CO/jhanwar_1982__moments_CO.csv", []).append(point)
        for marks, values, _text in _tokens(cells[2]):
            if len(values) == 2:  # Kirkwood-Unsold pair, footnote e (docstring)
                for stem, value in zip(
                    ("mulder_1980_unsold", "mulder_1980_kirkwood"), values, strict=True
                ):
                    path = f"CO/jhanwar_1982__{stem}_CO.csv"
                    files.setdefault(path, []).append(
                        Point(order, order, value, "", printed)
                    )
                continue
            for mark in marks.split(","):
                path = f"CO/jhanwar_1982__{CO_MARKS[mark]}_CO.csv"
                point = Point(order, order, values[0], "", printed)
                files.setdefault(path, []).append(point)
    return files


# Table 4's literature cells: how many entries are on the first (UV) line of each
# row; the rest are on the second (total) line.
CO2_FIRST_LINE = {"2": 1, "1": 1, "-1": 1, "-2": 4, "-4": 1, "-6": 1, "-8": 1}


def _co2_sums() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = table(DIRECTORY, TABLE_4)[1:]
    for order, cells in zip(ORDERS, rows, strict=True):
        printed = " | ".join(cells)
        # Recommended: "UV" then "total", joined by the extraction without a
        # separator ("14.1414.24"); every mantissa in the column has two decimals
        # (image), which splits them.
        numbers = [found[0] for found in _RECOMMENDED.finditer(cells[1])]
        assert len(numbers) == 2, cells[1]
        for part, text in zip(("uv", "total"), numbers, strict=True):
            path = f"CO2/jhanwar_1982__moments_{part}_CO2.csv"
            point = Point(order, order, _value(text), "", printed)
            files.setdefault(path, []).append(point)
        entries = _tokens(cells[2])
        first = CO2_FIRST_LINE.get(order, 0)
        for index, (marks, values, _text) in enumerate(entries):
            part = "uv" if index < first else "total"
            if len(values) == 2:  # Kirkwood-Unsold pair (docstring)
                for stem, value in zip(
                    ("mulder_1980_unsold", "mulder_1980_kirkwood"), values, strict=True
                ):
                    path = f"CO2/jhanwar_1982__{stem}_CO2.csv"
                    files.setdefault(path, []).append(
                        Point(order, order, value, "", printed)
                    )
                continue
            stem = CO2_MARKS[marks]
            suffix = f"_{part}" if stem in {"pack_1974", "watson_1936"} else ""
            path = f"CO2/jhanwar_1982__{stem}{suffix}_CO2.csv"
            point = Point(order, order, values[0], "", printed)
            files.setdefault(path, []).append(point)
    return files


def _check_pseudospectra() -> None:
    """Assert that the record's Table 5 terms are the table's cells."""
    record = tomllib.loads(RECORD.read_text(encoding="utf-8"))
    formulas = {f["id"]: f for f in record["formulas"]}
    rows = json_table(DIRECTORY, TABLE_5)[2:]
    for gas, (e_column, f_column) in {"CO": (0, 1), "CO2": (2, 3)}.items():
        terms = formulas[f"pseudospectrum_{gas}"]["terms"]
        assert terms["E"] == [_value(r[e_column]) for r in rows], gas
        assert terms["f"] == [_value(r[f_column]) for r in rows], gas


def extract() -> dict[str, list[Point]]:
    """Return the DOSDs (Tables 1-2) and the dipole sums (Tables 3-4)."""
    _check_pseudospectra()
    files = _dosd(TABLE_1, {2: "CO/jhanwar_1982__dosd_CO.csv"})
    files |= _dosd(
        TABLE_2,
        {
            2: "CO2/jhanwar_1982__dosd_uv_CO2.csv",
            3: "CO2/jhanwar_1982__dosd_total_CO2.csv",
        },
    )
    return files | _co_sums() | _co2_sums()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
