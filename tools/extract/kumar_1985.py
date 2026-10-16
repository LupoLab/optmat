# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar & Meath (1985): rare-gas integrated oscillator strengths, dipole sums.

The tables print numbers as a mantissa followed by the power of ten in parentheses,
``2.875(-1)`` = 2.875e-1, ``1.108(1)`` = 11.08, and (Tables 8-11) the reference of the
oscillator-strength data after it, ``1.528(-3)(66)`` or ``1.050 (70)``; the reference is
a two- or three-digit number or ``STOBBE``, the exponent a single signed digit written
directly after the mantissa. Values are written in e-notation (``2.875e-1``), the
printed cell kept in the ``printed`` column. The hydrogen halides (Tables 1-7 and the
HF/HCl/HBr columns) are not curated species and are skipped.

- Tables 8-11 (Ne ``tbl_9594b1549318``, Ar ``tbl_be26a24107d8``, Kr
  ``tbl_d6de19f0cf8e``, Xe ``tbl_b6cd2000760f``; cells files; images checked): the
  integrated dipole oscillator strength of the adopted DOSD over each energy interval
  (eV), first row the whole discrete spectrum (range in parentheses). Only the
  "Adopted DOSD" column is extracted; the columns of initial input data, one
  reference per cell, are not. The format has no quantity for an oscillator strength
  integrated over an interval: the files use ``oscillator_strength`` with an index
  abscissa (the row's position in the table) and the interval, in eV, as the row id
  (``21.47-21.57``, ``discrete:16.67-21.48``, ``100000-inf``). The adopted column of
  each table sums to the number of electrons (10.0003, 17.9994, 35.9987, 53.9976).
- Table 12 (``tbl_1b74f66df9be``, cells file; image checked): the dipole sums S_k,
  k = 2, 1, 0, -1/2, ..., -10, "Adopted" and "Leonard and Barker" per gas, rows by k
  (``-0.5`` for -1/2); '—' cells skipped. The logarithmic sums L_k are not extracted
  (no quantity).
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Kumar_1985_IntegratedDipoleOscillator"
DOSD_TABLES = {
    "Ne": "tbl_9594b1549318",
    "Ar": "tbl_be26a24107d8",
    "Kr": "tbl_d6de19f0cf8e",
    "Xe": "tbl_b6cd2000760f",
}
TABLE_12 = "tbl_1b74f66df9be"
# Table 12 column -> (gas, dataset stem).
TABLE_12_COLUMNS = {
    1: ("Ne", "moments"),
    2: ("Ne", "leonard_barker"),
    3: ("Ar", "moments"),
    4: ("Ar", "leonard_barker"),
    5: ("Kr", "moments"),
    6: ("Kr", "leonard_barker"),
    7: ("Xe", "moments"),
    8: ("Xe", "leonard_barker"),
}
_CELL = re.compile(
    r"^\s*(?P<mantissa>-?\d+\.\d+)\s*(?:\(\s*(?P<exponent>-?\d)\))?"
    r"\s*(?:\((?P<reference>\d{2,3}|STOBBE)\))?\s*$"
)
_ORDER = re.compile(r"^S_\{(?P<k>-?\d+(?:/2)?)\}")
_DISCRETE = re.compile(r"\(([\d.]+)-([\d.]+)\)")


def _value(cell: str) -> str:
    """Return a ``mantissa(exponent)(reference)`` cell's number in e-notation."""
    match = _CELL.match(cell.replace("\N{MINUS SIGN}", "-"))
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = match["mantissa"], match["exponent"]
    return f"{mantissa}e{exponent}" if exponent else mantissa


def _bound(text: str) -> str:
    """Normalise one end of a printed energy interval (eV)."""
    text = text.replace(" ", "").replace("\N{INFINITY}", "inf")
    return {"\\infty": "inf", "10^5": "100000"}.get(text, text)


def _interval(label: str) -> str | None:
    """Return the row id for an energy-interval label, or None for other rows."""
    found = _DISCRETE.search(label)
    if found:
        return f"discrete:{found[1]}-{found[2]}"
    if "-" not in label:
        return None  # the header and the bare "Discrete" label row
    low, high = label.split("-", 1)
    return f"{_bound(low)}-{_bound(high)}"


def _dosd(table_id: str) -> list[Point]:
    points = []
    for cells in json_table(DIRECTORY, table_id)[1:]:
        row = _interval(cells[0].strip())
        if row is None or not cells[1].strip():
            continue
        low, high = row.removeprefix("discrete:").split("-")
        printed = f"{cells[0]} | {cells[1]}"
        points.append(
            Point(row, low, _value(cells[1]), "", printed, extra={"x_upper": high})
        )
    return points


def _order(label: str) -> str | None:
    found = _ORDER.match(label)
    if found is None:
        return None  # header rows and the L_k rows
    k = found["k"]
    if k.endswith("/2"):
        return str(int(k.removesuffix("/2")) / 2)
    return k


def _moments() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_12):
        order = _order(cells[0].strip())
        if order is None:
            continue
        printed = " | ".join(cells)
        for column, (gas, stem) in TABLE_12_COLUMNS.items():
            if cells[column].strip() in {"\N{EM DASH}", ""}:
                continue
            point = Point(order, order, _value(cells[column]), "", printed)
            path = f"{gas}/kumar_1985__{stem}_{gas}.csv"
            files.setdefault(path, []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the adopted integrated DOSDs and the dipole sums of Table 12."""
    files = {
        f"{gas}/kumar_1985__dosd_{gas}.csv": _dosd(table_id)
        for gas, table_id in DOSD_TABLES.items()
    }
    files |= _moments()
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
