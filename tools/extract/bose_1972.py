# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bose, Sochanski & Cole (1972): dielectric virial coefficients of CH4.

Table I (``tbl_36d8a70060d0``, cells file only; image checked, the cells agree) gives
A_eps (cm^3 mol^-1), B_eps (cm^6 mol^-2) and C_eps (cm^9 mol^-3) of CH4 and CF4 at
279.8, 322.5 and 373.4 K, each printed ``value±error`` (three standard deviations,
p. 3). Only CH4 is extracted (CF4 is not a curated species). The gas cell spans three
rows in the image (merged cell), so the CH4 rows are the three after the header.

Structural fix: the C_eps header extracts as "cm9 mol-8"; the table image prints
cm^9 mol^-3, the unit of a third virial coefficient (the record's unit is cm9/mol3).

A_eps goes to one file with the temperature (K) as the row id; B_eps and C_eps to one
file per temperature (one temperature per dataset).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bose_1972_DielectricPressureVirial"
TABLE_I = "tbl_36d8a70060d0"
_PM = "\N{PLUS-MINUS SIGN}"
_HEADER = [
    "Gas",
    "T(°K)",
    r"A_{\epsilon}(cm3 mol-1)",
    r"B_{\epsilon}(cm6 mol-2)",
    # The image prints cm^9 mol^-3 (see the module docstring).
    r"C_{\epsilon}(cm9 mol-8)",
]


def split(cell: str) -> tuple[str, str]:
    """Return the value and the error of ``6.538±0.003``, as printed."""
    value, sign, sigma = cell.partition(_PM)
    if not sign:
        msg = f"no error in {cell!r}"
        raise ValueError(msg)
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return A_eps (one file) and B_eps, C_eps (one file per temperature) of CH4."""
    rows = json_table(DIRECTORY, TABLE_I)
    if rows[0] != _HEADER:
        msg = f"unexpected header {rows[0]!r}"
        raise ValueError(msg)
    if rows[1][0] != "CH_4" or rows[4][0] != "CF_4":
        msg = "the CH4 block is not rows 1-3"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    a_eps = []
    for cells in rows[1:4]:
        temperature = number(cells[1])
        printed = " | ".join(cells)
        a_eps.append(Point(temperature, "inf", *split(cells[2]), printed))
        files[f"CH4/bose_1972__B_eps_{temperature}_CH4.csv"] = [
            Point(temperature, "inf", *split(cells[3]), printed)
        ]
        files[f"CH4/bose_1972__C_eps_{temperature}_CH4.csv"] = [
            Point(temperature, "inf", *split(cells[4]), printed)
        ]
    files["CH4/bose_1972__A_eps_CH4.csv"] = a_eps
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
