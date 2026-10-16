# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Egan & Stone (2011): refractivity of N2 and Ar at 100 kPa, 20 C, 633 nm.

Table 2 (``tbl_7d11225de36b``, cells file only; image checked) gives (n - 1) x 10^8 for
nitrogen and argon at p = 100 kPa and T = 20 C: this work, and the values the authors
derive from the molar refractivities of Birch, Hou and Thalmann, and Achtermann et al.
One file per row and gas; the "+-" uncertainty goes to the sigma column.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Egan_2011_AbsoluteRefractometryDry"
TABLE = "tbl_7d11225de36b"
GASES = {1: "N2", 2: "Ar"}
ROWS = {
    "This work": "refractivity",
    "Birch": "birch_1991",
    "Hou and Thalmann": "hou_1994",
    "Achtermann et al.": "achtermann_1991",
}
WAVELENGTH = "632.9908"  # the He-Ne line; see the source record


def extract() -> dict[str, list[Point]]:
    """Return one file per row of Table 2 and gas."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE)[2:]:
        name = ROWS[cells[0].strip()]
        printed = " | ".join(cells)
        for column, gas in GASES.items():
            value, error = cells[column].split("\\pm")
            path = f"{gas}/egan_2011__{name}_{gas}.csv"
            files[path] = [
                Point(WAVELENGTH, WAVELENGTH, number(value), number(error), printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
