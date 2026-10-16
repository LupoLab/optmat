# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Smith, Parkinson & Huber (1975): Kr n - 1 at 168-288 nm (Table 3).

Table 3 has no CSV rendering (its two-level heading spans columns) and is read from its
JSON cells; the cells agree with the table image. Columns: vacuum wavelength (nm), this
work, the three-term Sellmeier interpolation ("calc."), their difference, Abjean et
al.'s data and difference, and Kronjager's and Koch's data as adjusted by Leonard and
difference. The difference columns are not extracted. All values ``(n - 1) x 10^6``.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Smith_1975_RefractiveIndexKrypton"
TABLE_3 = "tbl_bb31b7ebfefd"
COLUMNS = {
    1: "experiment_Kr",
    2: "interpolated_Kr",
    4: "abjean_1971_Kr",
    6: "kronjager_koch_Kr",
}


def extract() -> dict[str, list[Point]]:
    """Return one file per value column of Table 3."""
    rows = json_table(DIRECTORY, TABLE_3)[2:]  # two heading rows
    files: dict[str, list[Point]] = {}
    for column, stem in COLUMNS.items():
        files[f"Kr/smith_1975__{stem}.csv"] = [
            Point(cells[0], cells[0], number(cells[column]), "", " | ".join(cells))
            for cells in rows
            if cells[column].strip()
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
