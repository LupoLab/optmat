# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kumar & Thakkar (2011), supplement: integrated DOS of O2 and O3.

Page 2 (``tbl_528da2c47d04``, cells file only; image checked, every cell agrees): the
final O2 DOSD, columns j, E_j (eV), IDOS over [E_j, E_{j+1}] (the last interval ends at
infinity) and the source of the initial data; the first row, "Discrete", is the summed
strength of the discrete Schumann-Runge bands (Bethke; Lewis et al.). Rows j = 1-3 are
not printed (the discrete row stands for them).

Page 3 (``tbl_4d1cc5fd6d9e``, CSV; image checked, every cell agrees): the five O3 DOSDs
A, B1, B2, C and D, columns j, E(j) (eV, as for O2) and one IDOS column per
distribution, 39 intervals. The sources of the initial data per j are listed on page 4
as text (in the source record; one line, "j=14-15: Mason et al.", is missing from the
library's Markdown and survives only as an unclassified block of ``document.json``; the
page image shows it).

"The numbers in parentheses are powers of ten": ``8.132(-8)`` is written 8.132e-8, the
printed cell kept. The abscissa of an interval is its lower edge E_j (eV), as in
``kumar_2010_supplement``. The pseudospectra (page 5) are recorded as formulas in the
source record, not extracted here.
"""

import dataclasses
import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Kumar_2011_OzoneUnresolvedDiscrepancies/supplement_01"
TABLE_O2 = "tbl_528da2c47d04"
TABLE_O3 = "tbl_4d1cc5fd6d9e"
DISTRIBUTIONS = ("A", "B1", "B2", "C", "D")
_CELL = re.compile(r"^(?P<mantissa>\d+\.\d+)(?:\((?P<exponent>-?\d+)\))?$")


def _value(cell: str) -> str:
    """Return a ``mantissa(power of ten)`` cell in e-notation."""
    match = _CELL.match(cell.strip())
    if match is None:
        msg = f"unreadable cell {cell!r}"
        raise ValueError(msg)
    mantissa, exponent = match["mantissa"], match["exponent"]
    return f"{mantissa}e{exponent}" if exponent else mantissa


def _o2() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_O2)
    if rows[0][:3] != ["j", "Ej(eV)", "IDOS"]:
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    prefix = "O2/kumar_2011_supplement__"
    discrete, idos = [], []
    for cells in rows[1:]:
        label, energy, value, source = (c.strip() for c in cells)
        printed = " | ".join((label, energy, value, source))
        if label == "Discrete":
            discrete.append(Point("discrete", "1", _value(value), "", printed))
        else:
            idos.append(Point(label, _value(energy), _value(value), "", printed))
    return {
        f"{prefix}discrete_O2.csv": discrete,
        f"{prefix}idos_O2.csv": _with_upper_edges(idos),
    }


def _o3() -> dict[str, list[Point]]:
    header, *rows = table(DIRECTORY, TABLE_O3)
    if header != ["j", "E(j)", *DISTRIBUTIONS]:
        msg = f"unexpected header {header!r}"
        raise ValueError(msg)
    prefix = "O3/kumar_2011_supplement__idos_"
    files: dict[str, list[Point]] = {f"{prefix}{d}.csv": [] for d in DISTRIBUTIONS}
    for cells in rows:
        label, energy, *values = (c.strip() for c in cells)
        printed = " | ".join(c.strip() for c in cells)
        for name, value in zip(DISTRIBUTIONS, values, strict=True):
            files[f"{prefix}{name}.csv"].append(
                Point(label, _value(energy), _value(value), "", printed)
            )
    return files


def _with_upper_edges(points: list[Point]) -> list[Point]:
    """Close each interval at the next row's lower edge E_(j+1), the last at infinity.

    E_j is the lower edge of interval j, as in kumar_2010_supplement (checked there
    against kumar_1985's explicit intervals).
    """
    edges = [p.x for p in points[1:]] + ["inf"]
    return [
        dataclasses.replace(p, extra={"x_upper": edge})
        for p, edge in zip(points, edges, strict=True)
    ]


def extract() -> dict[str, list[Point]]:
    """Return the O2 and O3 IDOS intervals and the O2 discrete sum."""
    o3 = {path: _with_upper_edges(points) for path, points in _o3().items()}
    return _o2() | o3


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
