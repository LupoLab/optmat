# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Hunt, Poll and Wolniewicz (1984): polarizability matrix elements of H2, D2.

Tables 5-10 (p. 1721-1723; no CSVs, the JSON cells read) give the matrix elements
<0J|alpha_par|vJ'> and <0J|alpha_perp|vJ'> in atomic units between the ground level
(v = 0, J) and the levels v = 0-5, J' = J (Tables 5-6, J = 0-5), J + 2 (Tables 7-8,
J = 0-3) and J - 2 (Tables 9-10, J = 2-5), for H2, HD, HT, D2, DT and T2. Only the H2
and D2 columns are written, one file per table, gas and final v, one row per J.

Every cell of the H2 and D2 columns was compared with the table images (2026-09-27).
Structural points:

- The v label is a cell spanning its block's rows (empty below its first row in the
  cells); the J label is the second column. The D2 column of Table 5 is headed ``v_2``
  in the cells (``D_2`` on the page).
- Table 7's v = 1 block: the page prints the block label as ``1 0`` (v = 1, J = 0) with
  a stray mark, read as ``1-0`` spanning the v and J columns; in its rows J = 1-3 the J
  label moved into the H2 column and every value one column to the right, losing the T2
  value. On the page the rows read J = 1: 0.92915 (H2) ... 0.82238 (D2); J = 2: 0.84652
  ... 0.77311; J = 3: 0.76646 ... 0.72488, which is what the shifted cells hold.
- Table 5 prints D2's <00|alpha_par|20> as 0.05292 without a minus sign (image
  checked): a misprint, corrected by an erratum in the record, not here.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Hunt_1984_AbInitioCalculation"
# Table id -> (table number, polarizability component, Delta J = J' - J).
TABLES = {
    "tbl_d342fdc7648a": (5, "par", "dJ0"),
    "tbl_67070dd33681": (6, "perp", "dJ0"),
    "tbl_96e803fab9c9": (7, "par", "dJp2"),
    "tbl_9fd840e5f3aa": (8, "perp", "dJp2"),
    "tbl_69944584e92c": (9, "par", "dJm2"),
    "tbl_62b45dba2e04": (10, "perp", "dJm2"),
}
HEADERS = {
    5: ("V", "J", "H_2", "HD", "HT", "v_2", "DT", "T_2"),
}
HEADER = ("V", "J", "H_2", "HD", "HT", "D_2", "DT", "T_2")
COLUMNS = {"H2": 2, "D2": 5}
SHIFTED = "tbl_96e803fab9c9"  # Table 7, whose v = 1 block is shifted (docstring)


def _rows(table_id: str, number_: int) -> list[tuple[str, str, list[str]]]:
    """Return (v, J, cells) for every data row, the cells aligned with the header."""
    grid = json_table(DIRECTORY, table_id)
    if tuple(grid[0]) != HEADERS.get(number_, HEADER):
        msg = f"Table {number_} header {grid[0]!r}"
        raise ValueError(msg)
    rows = []
    v = ""
    for cells in grid[1:]:
        if table_id == SHIFTED and cells[0] == "1-0":
            # The v = 1, J = 0 row: 'v J' printed as one label.
            v, j, aligned = "1", "0", cells
        elif table_id == SHIFTED and v == "1" and not cells[0] and not cells[1]:
            # A shifted row of the same block: J in the H2 column, values one right.
            j, aligned = cells[2], ["", cells[2], *cells[3:], ""]
        else:
            v = cells[0] or v
            j, aligned = cells[1], cells
        if not (v.isdigit() and j.isdigit()):
            msg = f"Table {number_}: row {cells!r} has no v, J"
            raise ValueError(msg)
        rows.append((v, j, aligned))
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 matrix elements of Tables 5-10."""
    files: dict[str, list[Point]] = {}
    for table_id, (number_, component, delta) in TABLES.items():
        for v, j, cells in _rows(table_id, number_):
            printed = f"Table {number_} | v {v} | J {j} | " + " | ".join(cells[2:])
            for gas, column in COLUMNS.items():
                name = f"{gas}/hunt_1984__{component}_{delta}_v{v}_{gas}.csv"
                row = f"J={j}"
                files.setdefault(name, []).append(
                    Point(row, row, number(cells[column]), "", printed)
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
