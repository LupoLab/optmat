# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Lamy (1977): n and k of crystalline and fused quartz, 1000-1650 A (Table I).

Table I (``tbl_1c325c8f2d92``, cells only, no CSV) gives, per wavelength in angstrom,
N and K of crystalline quartz (ordinary ray) and of fused quartz, from a
Kramers-Kronig analysis of near-normal-incidence reflectance. Four data files: n and k
of each material. Blank K cells (crystalline quartz from 1450 A, fused quartz from 1550
A) are not printed values and are skipped.

Checked against the table image (2026-09-27), every cell. Two things the image shows:

- the crystalline-quartz K at 1200 A is printed "1 25": the decimal point is missing
  from the print (a gap where the point should be; the column's other cells, 1.45
  above and 1.05 below, print two decimals). The extracted cell has the same gap.
  Read as 1.25 (a structural reading of a lost point, not a changed digit); the
  ``printed`` column keeps "1 25";
- four small K values are printed with a third digit in parentheses, "0.08(4)",
  "0.05(5)", "0.08(0)", "0.04(3)": the number is read with that digit (0.084 ...);
  the paper does not say what the parentheses mean (presumably a less certain digit).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Lamy_1977_OpticalConstantsCrystalline"
TABLE = "tbl_1c325c8f2d92"
_PARENTHESISED = re.compile(r"^(\d+\.\d+)\((\d)\)$")


def _value(cell: str) -> str:
    """Return a cell's number: the lost point of "1 25", a parenthesised last digit."""
    if cell == "1 25":
        # Structural: the print lost the decimal point (see the docstring).
        return "1.25"
    match = _PARENTHESISED.match(cell)
    if match:
        return match.group(1) + match.group(2)
    return number(cell)


def extract() -> dict[str, list[Point]]:
    """Return n and k of crystalline (ordinary ray) and fused quartz."""
    rows = json_table(DIRECTORY, TABLE)
    if rows[0][0] != "Wavelength (Å)" or rows[1][1:] != ["N", "K", "N", "K"]:
        msg = f"unexpected header {rows[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {
        "SiO2_quartz/lamy_1977__quartz_o_n.csv": [],
        "SiO2_quartz/lamy_1977__quartz_o_k.csv": [],
        "SiO2/lamy_1977__fused_n.csv": [],
        "SiO2/lamy_1977__fused_k.csv": [],
    }
    names = list(files)
    for cells in rows[2:]:
        wavelength = number(cells[0])
        printed = " | ".join(cells)
        for name, cell in zip(names, cells[1:], strict=True):
            if cell:
                files[name].append(
                    Point(wavelength, wavelength, _value(cell), "", printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
