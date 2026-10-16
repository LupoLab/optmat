# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rice, Scuseria, Lee, Taylor & Almlof (neon manuscript): static gamma of Ne.

Table 1 (``tbl_b31243b9ed05``, cells file only; image checked, every cell agrees):
alpha and gamma of neon by basis set and method. The basis-set cell is printed once per
block (blank below it): it is carried down. Only gamma is taken. Row ids are "<basis>
<method>" with the footnote marks stripped; the two uncontracted bases differ only in
their footnote (b: 8 electrons correlated, c: 10 electrons correlated), so their rows
get " [8e]" and " [10e]".
"""

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Rice_1992_ConnectedTripleExcitations"  # formerly Unverified_bb10c94322e6
TABLE = "tbl_b31243b9ed05"
FILE = "Ne/rice_1992_neon__table1.csv"


def extract() -> dict[str, list[Point]]:
    """Return Table 1's gamma column."""
    grid = json_table(DIRECTORY, TABLE)
    if grid[0] != [
        "Basis set",
        "Method",
        "\N{GREEK SMALL LETTER ALPHA}",
        "\N{GREEK SMALL LETTER GAMMA}",
    ]:
        msg = f"unexpected header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    basis = ""
    for cells in grid[1:]:
        if cells[0]:
            basis = cells[0]
        name = basis.replace("^a", "")
        if name.endswith("^b"):
            name = name.removesuffix("^b") + " [8e]"
        elif name.endswith("^c"):
            name = name.removesuffix("^c") + " [10e]"
        name = name.replace(" ", "")
        value = cells[3]
        float(value)
        points.append(Point(f"{name} {cells[1]}", "inf", value, "", " | ".join(cells)))
    rows = [p.row for p in points]
    if len(set(rows)) != len(rows):
        msg = f"duplicate rows {rows!r}"
        raise ValueError(msg)
    return {FILE: points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
