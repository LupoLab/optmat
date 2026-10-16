# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Kalfoglou & Miller (1967): the He-Ar interaction second virial coefficient.

Table II (``tbl_4e0de8840ac9``, p. 1258, cells file; image checked) prints, for the
present work, B_He-Ar at 30, 100, 200, 300, 400 and 500 C, in atm^-1 x 10^4 and in
cc/mole, each "value +- probable error", and Tanner and Masson's values at 30 and 100 C.
The cc/mole column is extracted (the pair is filed under Ar, the partner He), rows by
the temperature in C; the atm^-1 column (a pressure-series coefficient, B/RT, whose
unit the vocabulary does not hold) and Tanner and Masson's values (another source's)
stay in ``printed``.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Kalfoglou_1967_CompressibilityGasesMixtures"
TABLE_II = "tbl_4e0de8840ac9"
TEMPERATURES = ("30", "100", "200", "300", "400", "500")


def split(cell: str) -> tuple[str, str]:
    """Return the value and the uncertainty of ``18.77 ± 0.10``."""
    value, separator, sigma = cell.partition("\N{PLUS-MINUS SIGN}")
    if not separator:
        msg = f"not a value +- uncertainty cell: {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return B_He-Ar in cc/mole against temperature."""
    grid = json_table(DIRECTORY, TABLE_II)
    if grid[1][:3] != ["Temp, °C", "Atm-1\N{MULTIPLICATION SIGN} 104", "cc/mole"]:
        msg = f"unexpected Table II header {grid[1]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[2:8]:
        temperature = number(cells[0])
        value, sigma = split(cells[2])
        points.append(
            Point(
                temperature,
                temperature,
                value,
                sigma,
                " | ".join(c for c in cells if c),
            )
        )
    if [p.row for p in points] != list(TEMPERATURES):
        msg = f"Table II temperatures {[p.row for p in points]!r}"
        raise ValueError(msg)
    return {"Ar/kalfoglou_1967__B12_He.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
