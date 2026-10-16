# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bideau-Mehu et al. (1976): n - 1 at 164.1 and 170.2 nm, revised krypton.

Both tables are uncaptioned in the library and have no CSV rendering (their headings
span columns); they are read from their JSON cells and were checked on the images.
``tbl_3c940398bc8b`` is Table 2: He, Ne, Ar, Xe and N2 at 164.1 and 170.2 nm, empty
cells being no measurement. ``tbl_d8d2bb41dadd`` is Table 3, krypton: the new values
and, reproduced, the authors' withdrawn values (Abjean et al. 1971), Smith et al.
(1975), Kronjager and Koch adjusted by Leonard (1974) and Chashchina and Shreider.
Both tables print ``(n - 1) x 10^6``.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "BideauMehu_1976_MeasurementsRefractiveIndexes"
TABLE_2 = "tbl_3c940398bc8b"
TABLE_3 = "tbl_d8d2bb41dadd"
TABLE_2_GASES = {1: "He", 2: "Ne", 3: "Ar", 4: "Xe", 5: "N2"}
TABLE_3_COLUMNS = {
    1: "experiment_Kr",
    2: "abjean_1971_Kr",
    3: "smith_1975_Kr",
    4: "leonard_1974_Kr",
    5: "chashchina_shreider_1969_Kr",
}


def _column(rows: list[list[str]], column: int) -> list[Point]:
    """Return a column's non-empty cells, the wavelength (nm) as row id."""
    return [
        Point(cells[0], cells[0], number(cells[column]), "", " | ".join(cells))
        for cells in rows
        if cells[column].strip()
    ]


def extract() -> dict[str, list[Point]]:
    """Return one file per gas of Table 2 and one per column of Table 3."""
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE_2)[2:]  # the factor row, the column heads
    for column, gas in TABLE_2_GASES.items():
        files[f"{gas}/bideau_mehu_1976__experiment_{gas}.csv"] = _column(rows, column)
    rows = json_table(DIRECTORY, TABLE_3)[2:]
    for column, stem in TABLE_3_COLUMNS.items():
        files[f"Kr/bideau_mehu_1976__{stem}.csv"] = _column(rows, column)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
