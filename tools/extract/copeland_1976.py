# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Copeland & Cole (1976): argon A_eps and the Bose & Cole reprints (Table I).

Table I (``tbl_a066d6e5d380``, CSV; image checked) lists the dielectric and pressure
virial coefficients of CF3H-Ar mixtures at 50 degrees C by the mole fraction of CF3H.
Only the pure-argon row (mole fraction 0) is taken. Its image prints two lines per
cell: this work, A_eps = 4.1406 ± 0.0014, and in parentheses the values of T. K. Bose
and R. H. Cole, J. Chem. Phys. 52, 140 (1970) (footnote a), A_eps = (4.1396 ±
0.0006) and B_eps = (0.72 ± 0.12), this work's B_eps being "...". The CSV joins the
two lines of a cell ("4.1406 ± 0.0014(4.1396 ± 0.0006)a", "... (0.72 ± 0.12)a");
they are split at the parenthesis. CF3H is not a curated species.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Copeland_1976_DielectricPressureVirial"
TABLE_I = "tbl_a066d6e5d380"
_PAIR = re.compile(r"^(?P<value>-?[\d.]+) ± (?P<error>[\d.]+)$")
_CELL = re.compile(r"^(?P<own>.*?)\s*\((?P<reprint>[^)]*)\)a$")


def pair(text: str) -> tuple[str, str]:
    """Return the value and the error of ``4.1406 ± 0.0014``, as printed."""
    match = _PAIR.match(text.strip())
    if match is None:
        msg = f"unreadable value {text!r}"
        raise ValueError(msg)
    return number(match["value"]), number(match["error"])


def extract() -> dict[str, list[Point]]:
    """Return argon's A_eps (this work) and Bose and Cole's A_eps and B_eps."""
    rows = table(DIRECTORY, TABLE_I)
    argon = rows[1]
    if argon[0] != "0":
        msg = f"row 1 is mole fraction {argon[0]!r}, not 0 (argon)"
        raise ValueError(msg)
    printed = " | ".join(argon)
    a_cell, b_cell = _CELL.match(argon[1]), _CELL.match(argon[2])
    if a_cell is None or b_cell is None or b_cell["own"] != "...":
        msg = f"unexpected argon cells {argon[1]!r}, {argon[2]!r}"
        raise ValueError(msg)
    return {
        "Ar/copeland_1976__A_eps_Ar.csv": [
            Point("323.15", "inf", *pair(a_cell["own"]), printed)
        ],
        "Ar/copeland_1976__A_eps_bose_cole_Ar.csv": [
            Point("Bose and Cole", "inf", *pair(a_cell["reprint"]), printed)
        ],
        "Ar/copeland_1976__B_eps_bose_cole_Ar.csv": [
            Point("Bose and Cole", "inf", *pair(b_cell["reprint"]), printed)
        ],
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
