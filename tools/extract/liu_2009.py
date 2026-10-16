# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Liu et al. (2009): the ionisation energy of H2 and earlier determinations.

Table III (``tbl_c65377ced290``, CSV; image checked, p. 8) prints E_i(H2), the
conventional (para-H2, N = 0 -> N+ = 0) ionisation energy in cm^-1, as determined over
forty years, one row per determination: the year, the value in the experiment or the
theory column, and the reference numbers of the paper's list (p. 9). The last row,
"2008 ... This work", is this paper's value; the others are other authors' values and
are written one file per row, as ``liu_2010`` does for D2.

Cells print digit-group spaces and the uncertainty in units of the last digit in
parentheses, "124 417.491 13(37)"; the spaces are dropped and the parenthesised digits
scaled to the last decimal. The ortho-H2 value exists only in the equation printed
under Table II (p. 6) and is an inline point of the record, as are Table II's
intervals (notes).
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Liu_2009_DeterminationIonizationDissociation"
TABLE_III = "tbl_c65377ced290"
_VALUE = re.compile(r"^([\d .]+?)(?:\((\d+)\))?$")
# (year, column, reference numbers as printed) -> file stem. The references are those
# of the paper's list (p. 9): 2 de Lange, Reinhold & Ubachs 2002; 4 Wolniewicz 1995;
# 24 Herzberg & Jungen 1972; 33 Herzberg 1969; 34 Eyler, Short & Pipkin 1986; 35 Glab &
# Hessler 1987; 36 Eyler et al. 1987; 37 McCormack et al. 1989; 38 Jungen et al. 1990;
# 39 Gilligan & Eyler 1992; 40 Kolos & Rychlewski 1993; 41 Shiner et al. 1993;
# 42 Wolniewicz 1993; 43 Meiners 1994 (M.S. thesis); 44 Kolos 1994.
ROWS = {
    ("1969", "expt", "33"): "herzberg_1969_expt",
    ("1969", "theory", "33"): "herzberg_1969_theory",
    ("1972", "expt", "24"): "herzberg_1972",
    ("1986", "expt", "34"): "eyler_1986",
    ("1987", "expt", "35"): "glab_1987",
    ("1987", "expt", "34, revised in 36"): "eyler_1987_ref34",
    ("1987", "expt", "35, revised in 36"): "eyler_1987_ref35",
    ("1989", "expt", "37"): "mccormack_1989",
    ("1990", "expt", "38"): "jungen_1990",
    ("1992", "expt", "39"): "gilligan_1992",
    ("1993", "theory", "40"): "kolos_1993",
    ("1993", "expt", "37, revised in 41"): "shiner_1993_ref37",
    ("1993", "expt", "38, revised in 41"): "shiner_1993_ref38",
    ("1993", "theory", "42"): "wolniewicz_1993",
    ("1994", "expt", "43"): "meiners_1994",
    ("1994", "theory", "44"): "kolos_1994",
    ("1995", "theory", "4"): "wolniewicz_1995",
    ("2002", "expt", "2"): "de_lange_2002",
    ("2008", "expt", "This work"): "ie_para",
}


def value_and_sigma(cell: str) -> tuple[str, str]:
    """Return the value and absolute uncertainty ("" if none) of a printed cell."""
    match = _VALUE.match(cell.strip())
    if match is None:
        msg = f"not a 'value(uncertainty)' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1).replace(" ", ""))
    if match.group(2) is None:
        return value, ""
    decimals = len(value.partition(".")[2])
    return value, f"{int(match.group(2)) * 10**-decimals:.{decimals}f}"


def extract() -> dict[str, list[Point]]:
    """Return the H2 ionisation energies of Table III, one file per row."""
    grid = table(DIRECTORY, TABLE_III)
    if grid[0] != ["Year", "Expt. (cm-1)", "Theor. (cm-1)", "Ref."]:
        msg = f"unexpected Table III header {grid[0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        year, expt, theory, reference = (cell.strip() for cell in cells)
        if bool(expt) == bool(theory):
            msg = f"Table III row {cells!r} is not one of experiment or theory"
            raise ValueError(msg)
        column = "expt" if expt else "theory"
        stem = ROWS[year, column, reference]
        value, sigma = value_and_sigma(expt or theory)
        files[f"H2/liu_2009__{stem}_H2.csv"] = [
            Point(year, "para-H2", value, sigma, " | ".join(cells))
        ]
    if len(files) != len(ROWS):
        msg = f"Table III rows missing: {len(ROWS) - len(files)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
