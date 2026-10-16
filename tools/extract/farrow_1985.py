# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Farrow & Rahn (1985): nonresonant CARS chi_nr of Ar and N2 (Table 2).

Table 2 (``tbl_0abf9a86f7ed``; the cells agree with the table image) compares chi_nr in
10^-18 cm^3/erg ('using the definition of relation (1)') for argon and nitrogen: this
work, Ref. 1 (Hall and Eckbreth's rescaling of Rado), Ref. 12 (Rosasco and Hurst),
Ref. 13 (Lundeen et al.) and Ref. 11 (Rado as rescaled by Lundeen et al.). Cells are
``value ± sigma`` or a bare value; "-" is no value.

Table 1 (``tbl_da579f77cd0b``, CSV, checked against the image): the argon-broadening
coefficients of the N2 Q(J) lines, HWHM x 10^3 cm^-1, printed in three (J, gamma_Ar)
column pairs as "37(3)" (the error in the last digit); J = 0 and 1 carry the mark "a"
("Extrapolated", the table's footnote, extracted as a page footnote). The J = 6 cell
reads "6'": the image shows a stray mark after the 6, not a prime; it is read as 6.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Farrow_1985_InterpretingCoherentAnti"
TABLE = "tbl_0abf9a86f7ed"
GASES = {"Argon": "Ar", "Nitrogen": "N2"}
# column -> (dataset prefix, pump wavelength in nm of the original measurement)
COLUMNS = {
    1: ("this_work", "532"),
    2: ("hall_eckbreth", "694.3"),
    3: ("rosasco_hurst", "676.4"),
    4: ("lundeen", "532"),
    5: ("rado_lundeen", "694.3"),
}


TABLE_1 = "tbl_da579f77cd0b"
_GAMMA = re.compile(r"^(?P<value>\d+)(?:\((?P<sigma>\d+)\))?(?:\^\{a\})?$")


def _table_1() -> list[Point]:
    """Return Table 1's argon-broadening coefficients, in J order."""
    grid = table(DIRECTORY, TABLE_1)
    if grid[0] != ["J", r"\gamma_{\text{Ar}}"] * 3:
        msg = f"unexpected Table 1 header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for pair in range(3):
        for cells in grid[1:]:
            j, cell = cells[2 * pair], cells[2 * pair + 1]
            if not j:
                continue
            match = _GAMMA.match(cell)
            if match is None:
                msg = f"unreadable Table 1 cell {cell!r}"
                raise ValueError(msg)
            row = f"Q({number(j.rstrip("'"))})"
            sigma = match["sigma"] or ""
            points.append(Point(row, row, match["value"], sigma, f"{j} | {cell}"))
    return points


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and Table 2 column, and the Table 1 file."""
    files = {"N2/farrow_1985__table1_Ar.csv": _table_1()}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, (prefix, wavelength) in COLUMNS.items():
            cell = cells[column].strip()
            if cell == "-":
                continue
            value, _, sigma = cell.partition("±")
            files[f"{gas}/farrow_1985__{prefix}_{gas}.csv"] = [
                Point(
                    "chi_nr",
                    wavelength,
                    number(value),
                    number(sigma) if sigma else "",
                    printed,
                )
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
