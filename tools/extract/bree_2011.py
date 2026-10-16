# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Bree, Demircan & Steinmeyer (2011): noble-gas n2 from a Kramers-Kronig model.

Table I (``tbl_3d9fe63acf83``) has no CSV (merged header): read from its JSON cells,
which equal the table image. Columns: n2 of the neutral atoms at 800 nm ("Z = 1"), of
the singly charged ions ("Z = 2", not curated: skipped), the static limit n2(w -> 0),
and "Refs. [18,21]", experimental data compiled from Bishop & Pipin (1989) and Shelton &
Rice (1994). All in 10^-8 cm^2/TW.

Per gas, the model's 800 nm and static values form one dataset (rows ``800nm`` and
``static``), the compiled reference column another (row ``refs_18_21``).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bree_2011_SaturationAllOptical"
TABLE = "tbl_3d9fe63acf83"
GASES = ("He", "Ne", "Ar", "Kr", "Xe")


def extract() -> dict[str, list[Point]]:
    """Return the model (800 nm, static) and compiled reference n2 per gas."""
    grid = json_table(DIRECTORY, TABLE)
    assert grid[1][1:] == [
        "Z = 1",
        "Z = 2",
        r"n_{2}(\omega \rightarrow 0)",
        "Refs. [18,21]",
    ]
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        gas = cells[0]
        assert gas in GASES, cells
        printed = " | ".join(cells)
        files[f"{gas}/bree_2011__kk_{gas}.csv"] = [
            Point("800nm", "800", number(cells[1]), "", printed),
            Point("static", "inf", number(cells[3]), "", printed),
        ]
        files[f"{gas}/bree_2011__refs_{gas}.csv"] = [
            Point("refs_18_21", "inf", number(cells[4]), "", printed)
        ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
