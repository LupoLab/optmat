# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Geints, Kompanets & Chekalin (2025): effective n2 of Ar, N2, CO2 at 800 nm.

Table 1 (``tbl_2c05f4da79d1``) has no CSV (merged header cells): read from its JSON
cells, which equal the table image. n2 in 10^-24 m^2/W at 1, 5 and 11 bar, each in the
short-pulse (SP, 50 fs) and long-pulse (LP, 250-500 fs chirped) regime; cells are
``9.7 ± 1.5^a`` with footnote a "Ref. [25]" (Wahlstrand et al. 2012) and b "This work";
"-" where Ar has no LP value. Row ids are ``<regime>_<pressure>bar``.

Files: the two reproduced 1 bar SP values (Ar, N2, footnote a) each on their own; this
work's values per gas and regime with the rows 1, 5 and 11 bar (the dataset excludes the
5 and 11 bar rows, see the source record).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Geints_2025_MeasurementsOpticalKerr"
TABLE = "tbl_2c05f4da79d1"
GASES = {"Ar": "Ar", "N_2": "N2", "CO_2": "CO2"}
COLUMNS = {
    1: ("sp", "1"),
    2: ("lp", "1"),
    3: ("sp", "5"),
    4: ("lp", "5"),
    5: ("sp", "11"),
    6: ("lp", "11"),
}
_CELL = re.compile(r"^([\d.]+) ± ([\d.]+)\^([ab])$")


def extract() -> dict[str, list[Point]]:
    """Return the reproduced and this work's values."""
    grid = json_table(DIRECTORY, TABLE)
    assert grid[0][1:] == ["1 bar", "", "5 bar", "", "11 bar", ""], grid[0]
    assert grid[1][1:] == ["SP, n_2^{sh}", "LP, n_2^{lng}"] * 3, grid[1]
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        for column, (regime, bar) in COLUMNS.items():
            if cells[column] == "-":
                continue
            match = _CELL.match(cells[column])
            assert match, cells[column]
            value, sigma, mark = match.groups()
            point = Point(
                f"{regime}_{bar}bar", "800", number(value), number(sigma), printed
            )
            if mark == "a":
                name = f"{gas}/geints_2025__wahlstrand_2012_{gas}.csv"
            else:
                name = f"{gas}/geints_2025__{regime}_{gas}.csv"
            files.setdefault(name, []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
