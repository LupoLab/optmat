# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar, Fairley & Meath (1985): the SF6 constrained DOSD and dipole sums.

The tables print powers of ten in parentheses, with spaces (``7.437( - 2)``,
``1.899 71( - 1)`` = 1.89971e-1); numbers are written in e-notation with the digits
unchanged, the cell kept in ``printed``. Images checked against the cells
(2026-09-27).

- Table I (``tbl_e5a10522c2c9``, cells file): the adopted (UV and higher energy) DOSD
  integrated over 27 energy regions (eV) from 5.0 eV to infinity, the "Adopted values"
  column only (its reference in the next cell is kept in ``printed``); the columns of
  initial data (one reference per cell) are not extracted. Row = the region, x and
  x_upper its edges (``10^5`` -> 100000, ``\infty`` -> inf).
- Table II (``tbl_bfe8065fdca1``, cells file): the adopted dipole sums, printed as
  two lines per property, the first "UV + higher energy" and the second "total
  including IR", in two column blocks (S_2 ... S_-5/2; S_-3, S_-4 and the mean
  excitation energies I_k, which have no quantity and are not extracted). Rows by k
  (``-0.5`` for -1/2, ...); the ``^a`` mark on the total S_-2 (a constraint) is
  stripped.
- Table III (``tbl_98bfede8c63d``): the pseudo-DOSDs are formulas in the record; this
  script asserts that the record's terms equal the table's cells.
"""

import re
import tomllib
from pathlib import Path

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Kumar_1985_DipolePropertiesDispersion"
TABLE_I = "tbl_e5a10522c2c9"
TABLE_II = "tbl_bfe8065fdca1"
TABLE_III = "tbl_98bfede8c63d"
RECORD = Path(__file__).resolve().parents[2] / (
    "src/optmat/data/sources/kumar_1985b.toml"
)
PREFIX = "SF6/kumar_1985b__"
_NUMBER = re.compile(r"(?P<mantissa>\d+\.\d+)(?:\((?P<exponent>-?\d+)\))?")
_ORDER = re.compile(r"^S_\{?(?P<k>-?\d+(?:/2)?)\}?$")
_EDGES = {"10^5": "100000", "\\infty": "inf"}


def _value(cell: str) -> str:
    """Return a cell's number in e-notation, the digits unchanged."""
    text = cell.replace(" ", "").replace("^a", "")
    found = _NUMBER.fullmatch(text)
    if found is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = found["mantissa"], found["exponent"]
    return number(f"{mantissa}e{exponent}" if exponent else mantissa)


def _order(label: str) -> str | None:
    found = _ORDER.match(label.strip())
    if found is None:
        return None
    k = found["k"]
    return str(int(k.removesuffix("/2")) / 2) if k.endswith("/2") else k


def _table_i() -> list[Point]:
    points = []
    for cells in json_table(DIRECTORY, TABLE_I)[1:]:
        low, high = (_EDGES.get(t, t) for t in cells[0].split("-", 1))
        row = f"{low}-{high}"
        printed = " | ".join(c for c in cells if c)
        points.append(Point(row, low, _value(cells[1]), "", printed, {"x_upper": high}))
    return points


def _table_ii() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_II)[1:]
    for block in (0, 3):  # the two column blocks: label, qualifier, value
        order = None
        for cells in rows:
            label = cells[block]
            if label:
                order = _order(label)
                part = "uv"
            else:
                part = "total"  # the property's second line (docstring)
            if order is None:
                continue  # a mean excitation energy I_k
            printed = " | ".join(cells[block : block + 3])
            point = Point(order, order, _value(cells[block + 2]), "", printed)
            files.setdefault(f"{PREFIX}moments_{part}_SF6.csv", []).append(point)
    return files


def _check_pseudospectra() -> None:
    """Assert that the record's Table III terms are the table's cells."""
    record = tomllib.loads(RECORD.read_text(encoding="utf-8"))
    formulas = {f["id"]: f for f in record["formulas"]}
    rows = json_table(DIRECTORY, TABLE_III)[2:]
    uv = formulas["pseudospectrum_SF6"]["terms"]
    assert uv["E"] == [_value(r[0]) for r in rows if r[0]]
    assert uv["f"] == [_value(r[1]) for r in rows if r[1]]
    total = formulas["pseudospectrum_SF6_total"]["coefficients"]
    for i, cells in enumerate(rows, start=1):
        assert total[f"E{i}"] == _value(cells[2]), i
        assert total[f"f{i}"] == _value(cells[3]), i


def extract() -> dict[str, list[Point]]:
    """Return the adopted DOSD (Table I) and the dipole sums (Table II)."""
    _check_pseudospectra()
    return {f"{PREFIX}dosd_SF6.csv": _table_i()} | _table_ii()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
