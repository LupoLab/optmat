# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Liu et al. (2010): the ionisation energies of ortho- and para-D2.

Table III (``tbl_fcc4c6682cca``, image checked, p. 7): the energy intervals (cm^-1)
whose sums (4) and (8) are the ionisation energies of ortho-D2 (N = 0 -> N+ = 0) and
para-D2 (N = 1 -> N+ = 1). Only those two rows are ionisation energies; the others
(and D0) are in the record's notes.

Table IV (``tbl_e28c58907bc7``, image checked, p. 8): E_i(D2) (ortho, N = 0 ->
N+ = 0) from the literature, one file per row (year and reference), experiment or
theory column. The 2010 row is this work's Table III value, and is not repeated.

Cells print digit-group spaces (in Table IV's first row a LaTeX thin space "\,") and
the uncertainty in units of the last digit in parentheses, "124 745.394 07(58)"; the
spaces are dropped and the parenthesised digits scaled to the last decimal. A trailing
footnote mark (Table IV's 1975 row, "^a": a -1.0 cm^-1 pressure-shift correction is
included) is removed from the value and kept in the ``printed`` column.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Liu_2010_DeterminationIonizationDissociation"
TABLE_3 = "tbl_fcc4c6682cca"
TABLE_4 = "tbl_e28c58907bc7"
_VALUE = re.compile(r"^([\d .]+?)(?:\((\d+)\))?(?:\^\{?a\}?)?$")
TABLE_3_ROWS = {
    "(4)=(1)+(2)+(3)": ("D2/liu_2010__ie_ortho_D2.csv", "ortho-D2"),
    "(8)=(5)+(6)+(7)": ("D2/liu_2010__ie_para_D2.csv", "para-D2"),
}
# (year, reference numbers as printed) -> file stem; the references are those of the
# paper's list (p. 12): 28 Takezawa & Tanaka 1975; 29 Jungen et al. 1990; 30, 31
# Kolos et al. 1986 and Wolniewicz & Orlikowski 1991; 32, 33 Gilligan & Eyler 1992 and
# Jungen et al. 1992; 34 Shiner et al. 1993; 35 Kolos & Rychlewski 1993; 36
# Wolniewicz 1993; 37 Kolos 1994; 38 Wolniewicz 1995.
TABLE_4_ROWS = {
    ("1975", "28"): "takezawa_1975",
    ("1990", "29"): "jungen_1990",
    ("1991", "30,31"): "kolos_1991",
    ("1992", "32,33"): "gilligan_1992",
    ("1993", "34"): "shiner_1993",
    ("1993", "35"): "kolos_1993",
    ("1993", "36"): "wolniewicz_1993",
    ("1994", "37"): "kolos_1994",
    ("1995", "38"): "wolniewicz_1995",
}


def value_and_sigma(cell: str) -> tuple[str, str]:
    """Return the value and absolute uncertainty ("" if none) of a printed cell."""
    text = cell.strip().replace("\\,", "")
    match = _VALUE.match(text)
    if match is None:
        msg = f"not a 'value(uncertainty)' cell: {cell!r}"
        raise ValueError(msg)
    value = number(match.group(1).replace(" ", ""))
    if match.group(2) is None:
        return value, ""
    decimals = len(value.partition(".")[2])
    return value, f"{int(match.group(2)) * 10**-decimals:.{decimals}f}"


def _table_3() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_3)[1:]:
        label = cells[0].strip()
        if label not in TABLE_3_ROWS:
            continue
        path, x = TABLE_3_ROWS[label]
        value, sigma = value_and_sigma(cells[2])
        files[path] = [Point(label, x, value, sigma, " | ".join(cells))]
    if len(files) != len(TABLE_3_ROWS):
        msg = f"Table III rows missing: {set(TABLE_3_ROWS) - set(files)}"
        raise ValueError(msg)
    return files


def _table_4() -> dict[str, list[Point]]:
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_4)[1:]:
        year, expt, theory, reference = (cell.strip() for cell in cells)
        if reference == "This work":
            continue  # Table III (4), extracted there
        stem = TABLE_4_ROWS[year, reference]
        if bool(expt) == bool(theory):
            msg = f"Table IV row {cells!r} is not one of experiment or theory"
            raise ValueError(msg)
        value, sigma = value_and_sigma(expt or theory)
        files[f"D2/liu_2010__{stem}_D2.csv"] = [
            Point(year, "ortho-D2", value, sigma, " | ".join(cells))
        ]
    if len(files) != len(TABLE_4_ROWS):
        msg = "Table IV rows missing"
        raise ValueError(msg)
    return files


def extract() -> dict[str, list[Point]]:
    """Return the D2 ionisation energies of Tables III and IV."""
    return _table_3() | _table_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
