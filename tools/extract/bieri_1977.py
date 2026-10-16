# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bieri et al. (1977): the valence ionisation energies of CH4.

The first table of Section A, Linear Hydrocarbons (``tbl_1342dbc497d6``, cells only;
image checked, p. 2218): CH4, C2H6, C2H4, C2H2, each with an orbital column, a value
column and a reference column. Only the CH4 columns are taken: the 1t2 band's three
Jahn-Teller maxima, printed in one cell as "13.6^{a)}14.4^{a)}15.0^{a)}", and the 2a1
band, "22.9^{b)}". The footnote marks (a: maxima of the Jahn-Teller distorted band; b:
the adiabatic energy of that band, in the footnote) are removed from the values and
kept in the ``printed`` column. The reference cell of the 1t2 row, "[17]", applies to
every value below it in the column (the legend, p. 2218), so all four are reproduced.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bieri_1977_ValenceIonizationEnergies"
TABLE = "tbl_1342dbc497d6"
_MAXIMUM = re.compile(r"([\d.]+)\^\{([ab])\)\}")
ORBITALS = {"1t_2": ("1t2", "a"), "2a_1": ("2a1", "b")}


def extract() -> dict[str, list[Point]]:
    """Return the band maxima of CH4."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[0][1] != "CH_4":
        msg = f"the first column group is not CH4: {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[2:]:
        orbital = cells[1].strip()
        if not orbital:
            continue
        label, mark = ORBITALS[orbital]
        printed = " | ".join(cells[:4])
        entries = _MAXIMUM.findall(cells[2])
        if not entries or any(found != mark for _, found in entries):
            msg = f"unexpected CH4 cell {cells[2]!r}"
            raise ValueError(msg)
        for index, (value, _) in enumerate(entries, start=1):
            row = f"{label}_{index}" if len(entries) > 1 else label
            points.append(Point(row, label, number(value), "", printed))
    return {"CH4/bieri_1977__vertical_CH4.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
