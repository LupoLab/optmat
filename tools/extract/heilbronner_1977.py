# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Heilbronner, Jones & Maier (1977): the sigma ionisation energies of CH4.

Table 2 (``tbl_d6b9168f4f46``, CSV; image checked, p. 1700) compares the sigma-type
ionisation energies (eV, "all values +-0.1 eV") of the alkanes RH with those of the
alkyl groups R in the diacetylenes R-(C#C)2-R, in column triples (alkane, R, the
difference). Only the first column, CH4, is a curated species: its three entries are
the maxima of CH4's Jahn-Teller split first band (the 1t2 ionisation; see the record).
The rows are labelled as ``bieri_1977`` labels the same three values (1t2_1 ...).
"""

from _library import Point, number, table, write

DIRECTORY = "Heilbronner_1977_IonizationEnergiesDi"
TABLE_2 = "tbl_d6b9168f4f46"


def extract() -> dict[str, list[Point]]:
    """Return CH4's three 1t2 band maxima from Table 2."""
    grid = table(DIRECTORY, TABLE_2)
    if grid[0][0] != "CH_4":
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    cells = [row[0].strip() for row in grid[1:] if row[0].strip()]
    if len(cells) != 3:
        msg = f"expected three CH4 values, found {cells!r}"
        raise ValueError(msg)
    points = [
        Point(f"1t2_{i}", "1t2", number(cell), "", f"CH_4 | {cell}")
        for i, cell in enumerate(cells, start=1)
    ]
    return {"CH4/heilbronner_1977__vertical_CH4.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
