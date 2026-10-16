# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Scrinzi, Geissler & Brabec (1999): the static-field ionisation rates of He.

Table I (``tbl_695d85562032``, CSV; image checked, p. 708) prints 60 pairs (E, w), the
field and the rate in atomic units, in three column pairs of 20 rows (E = 0.08-0.27,
0.28-0.47 and 0.48-1.70 au). They are written in field order as one file.

Structural fix: the extraction renders "0.463 \times 10^{-7}" as "0.463 3 1027": the
multiplication sign became "3" and the minus sign of the exponent "2" (the image shows
"0.463 x 10^-7"; every exponent in the table is negative, from -7 to -1). Such a cell is
rewritten "0.463e-7"; ``printed`` keeps the cell as extracted. Cells without a factor
(0.105 onwards) are plain numbers.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Scrinzi_1999_IonizationAboveCoulomb"
TABLE_I = "tbl_695d85562032"
_MANGLED = re.compile(r"^(\d\.\d+) 3 102(\d)$")


def rate(cell: str) -> str:
    """Return a printed rate as a plain number string, undoing the lost glyphs."""
    text = cell.strip()
    if match := _MANGLED.match(text):
        return f"{match.group(1)}e-{match.group(2)}"
    return number(text)


def extract() -> dict[str, list[Point]]:
    """Return He's 60 static rates, in field order."""
    grid = table(DIRECTORY, TABLE_I)
    if grid[0] != ["E", "w"] * 3:
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for pair in range(3):
        for cells in grid[1:]:
            field, printed = cells[2 * pair].strip(), cells[2 * pair + 1].strip()
            points.append(
                Point(
                    number(field),
                    number(field),
                    rate(printed),
                    "",
                    f"{field} | {printed}",
                )
            )
    if len(points) != 60 or [float(p.x) for p in points] != sorted(
        float(p.x) for p in points
    ):
        msg = "Table I fields are not 60 increasing values"
        raise ValueError(msg)
    return {"He/scrinzi_1999__table1_He.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
