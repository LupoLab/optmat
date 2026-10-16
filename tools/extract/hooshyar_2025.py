# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hooshyar & Dumitrache (2025): dissociative recombination coefficients.

Table 2 (``tbl_a903ef2bd6e1``, CSV; image checked, p. 8), 'Reprinted from [28]'
(Shneider, Zheltikov & Miles 2011): the air-plasma reactions of the model. Of the
dissociative recombinations, R12 (N2+), R13 (O2+) and R14 (NO+) concern registered
ions; each is printed "A x 10^-13 (300/T_e)^n" m^3 s^-1 (n = 1/2, 1, 1.5). The mantissa
is written (the record's scale holds 10^-13); the printed form, the exponent and the
unit are checked against the record's values here.

Register: TC-979.
"""

import re

from _library import Point, table, write

DIRECTORY = "Hooshyar_2025_ModelingFemtosecondFilament"
TABLE_2 = "tbl_a903ef2bd6e1"
# Row, the data file, and the exponent printed (as in the record's parameters).
ROWS = {
    "R12": ("N2+/hooshyar_2025__dr_N2+.csv", "^{1/2}"),
    "R13": ("O2+/hooshyar_2025__dr_O2+.csv", ""),
    "R14": ("NO+/hooshyar_2025__dr_NO+.csv", "^{1.5}"),
}
_RATE = re.compile(
    r"^([\d.]+) \\times 10\^\{-13\}"
    r"\\left\(\\frac\{300\}\{T_\{\\mathrm\{e\}\}\}\\right\)(.*)$"
)


def extract() -> dict[str, list[Point]]:
    """Return R12-R14's coefficients at 300 K."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_2):
        if cells[0] not in ROWS:
            continue
        path, exponent = ROWS[cells[0]]
        match = _RATE.match(cells[2])
        if match is None or match.group(2) != exponent or cells[3] != "[m^{3}s^{-1}]":
            msg = f"{cells[0]} printed as {cells[2]!r} {cells[3]!r}"
            raise ValueError(msg)
        files[path] = [Point(cells[0], "", match.group(1), "", " | ".join(cells))]
    if len(files) != len(ROWS):
        msg = f"found {sorted(files)}, not every row of {sorted(ROWS)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
