# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bischel and Dyer (1986, JOSA B): the absolute Raman gain of H2 Q(1).

Table 1 (``tbl_33f36ba412af``, p. 680; no CSV, the JSON cells match the image) compares
the measured steady-state gain coefficients (this work, and previous work with their
footnoted references) with Eq. (11), in 10^-9 cm/W at 298 K, per pump wavelength. A dash
marks an empty cell. Previous work is footnoted a (Ref. 8), b (Ref. 10), c (Ref. 12);
the 265 nm value has no footnote. Entries print ``value \pm error``.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bischel_1986_WavelengthDependenceAbsolute"
TABLE_1 = "tbl_33f36ba412af"
DASH = "—"
_ENTRY = re.compile(r"^([\d.]+) \\pm ([\d.]+)(\^\{[abc]\})?$")


def _entry(cell: str) -> tuple[str, str]:
    """Return the value and error of a ``value +- error`` cell (printed with pm)."""
    match = _ENTRY.match(cell.strip())
    if match is None:
        msg = f"unexpected entry {cell!r}"
        raise ValueError(msg)
    return number(match.group(1)), number(match.group(2))


def extract() -> dict[str, list[Point]]:
    """Return Table 1's three columns as data files."""
    measured: list[Point] = []
    previous: list[Point] = []
    calculated: list[Point] = []
    for cells in json_table(DIRECTORY, TABLE_1)[2:]:
        wavelength, printed = cells[0], " | ".join(cells)
        for cell, points in ((cells[1], measured), (cells[2], previous)):
            if cell != DASH:
                value, sigma = _entry(cell)
                points.append(Point(wavelength, wavelength, value, sigma, printed))
        calculated.append(Point(wavelength, wavelength, number(cells[3]), "", printed))
    return {
        "H2/bischel_1986b__gain_H2.csv": measured,
        "H2/bischel_1986b__gain_previous_H2.csv": previous,
        "H2/bischel_1986b__gain_eq11_H2.csv": calculated,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
