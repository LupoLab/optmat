# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Ditchburn & Heddle (1954): O2 Schumann-Runge band oscillator strengths.

Table 1 (``tbl_2b2a1247ae97``, CSV; image checked, the cells agree) prints the band
oscillator strengths f of the bands (v', 0), v' = 0-20, in three blocks of two rows (a
row "v'" of upper-state vibrational numbers, then a row "f" of values), each cell
"<mantissa> \times 10^{<k>}" with its own power of ten. The exponent varies from cell
to cell, so no column factor can be removed: the value column holds the number as
printed in E notation ("1.55e-8"), the printed column the raw cell. The text (p. 514)
says the values for v' = 1, 18, 19 and 20 are extrapolated (not marked in the table).

Row ids are the bands "(v',0)", the abscissa the index v'.
"""

import re

from _library import Point, table, write

DIRECTORY = "Ditchburn_1954_AbsorptionCrossSections"
TABLE_1 = "tbl_2b2a1247ae97"
PREFIX = "O2/ditchburn_1954__"

_F = re.compile(r"^(\d\.\d+) \\times 10\^\{(-\d+)\}$")


def _value(cell: str) -> str:
    match = _F.match(cell.strip())
    if match is None:
        msg = f"unexpected f cell {cell!r}"
        raise ValueError(msg)
    return f"{match.group(1)}e{match.group(2)}"


def extract() -> dict[str, list[Point]]:
    """Return the band oscillator strengths f(v', 0), v' = 0-20."""
    rows = table(DIRECTORY, TABLE_1)
    points = []
    for upper, strengths in zip(rows[0::2], rows[1::2], strict=True):
        if upper[0] != "v'" or strengths[0] != "f":
            msg = f"unexpected row labels {upper[0]!r}, {strengths[0]!r}"
            raise ValueError(msg)
        for v, f in zip(upper[1:], strengths[1:], strict=True):
            int(v)
            points.append(Point(f"({v},0)", v, _value(f), "", f"v' = {v} | f = {f}"))
    if [p.x for p in points] != [str(v) for v in range(21)]:
        msg = "the bands are not v' = 0-20 in order"
        raise ValueError(msg)
    return {f"{PREFIX}band_f.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
