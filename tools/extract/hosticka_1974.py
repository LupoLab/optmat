# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hosticka & Bose (1974): dielectric virial coefficients of SF6.

Table I (``tbl_23fb0fd23746``, no CSV: read from the JSON cells; image checked, the
cells agree) gives, for "This work" at 323.3, 348.3 and 373.9 K, A_eps (cm^3 mole^-1),
B_eps (cm^6 mole^-2), C_eps (cm^9 mole^-3) and B_p, each printed ``value ± error``
(three standard deviations, p. 1319), and "Literature" rows. A_eps goes to one file
(row = temperature in K); B_eps and C_eps to one file per temperature. The literature
row at 323.3 K (A_eps, B_eps, C_eps with the footnote mark "a", Reference 1: Nelson and
Cole 1971) is reproduced, one file per quantity; the other literature cells are
pressure virial coefficients B_p (no quantity; the record's notes).

The 323.3 K C_eps is printed "6900 ± 899" without the minus sign of the other two
temperatures and of the literature value; it is extracted as printed and corrected by
an erratum in the record.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hosticka_1974_DielectricPressureVirial"
TABLE_I = "tbl_23fb0fd23746"
_PM = "\N{PLUS-MINUS SIGN}"
_FOOTNOTE_A = re.compile(r"^(-?[\d.]+)a$")
QUANTITIES = {2: "A_eps", 3: "B_eps", 4: "C_eps"}


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the error of ``16.543 ± 0.017``."""
    value, sign, sigma = cell.partition(_PM)
    if not sign:
        msg = f"no error in {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return Table I's A_eps, B_eps and C_eps of SF6."""
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[0][:2] != ["Source", "T(K)"] or not grid[0][4].startswith("C_"):
        msg = f"unexpected Table I header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:7]:
        temperature = number(cells[1])
        printed = " | ".join(cells)
        if cells[0] == "This work":
            for column, name in QUANTITIES.items():
                point = Point(temperature, "inf", *_split(cells[column]), printed)
                stem = "A_eps_SF6" if name == "A_eps" else f"{name}_{temperature}_SF6"
                files.setdefault(f"SF6/hosticka_1974__{stem}.csv", []).append(point)
        elif cells[0] == "Literature" and cells[2]:
            for column, name in QUANTITIES.items():
                match = _FOOTNOTE_A.match(cells[column])
                if match is None:
                    msg = f"unexpected literature cell {cells[column]!r}"
                    raise ValueError(msg)
                point = Point(temperature, "inf", number(match[1]), "", printed)
                files[f"SF6/hosticka_1974__{name}_nelson_1971_SF6.csv"] = [point]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
