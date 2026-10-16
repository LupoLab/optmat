# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Christodouleas et al. (2009): static gamma components of ozone (Table 2).

Table 2 (``tbl_c26b470b9b60``, cells file only; image checked, every O3 cell agrees):
dipole moment, alpha, beta and gamma components of O3, SO2 and their ring isomers by
CCSD(T) and the worst and best DFT functional. Only the open (C2v) ozone rows and the
six gamma components are taken (one file per component; rows the methods). SO2 is not a
curated species; the ring (D3h) ozone is not the gas's ground state; alpha is linear
data.
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Christodouleas_2009_TrendsBondingEffect"
TABLE_2 = "tbl_c26b470b9b60"
COMPONENTS = {9: "xxxx", 10: "yyyy", 11: "zzzz", 12: "xxyy", 13: "yyzz", 14: "xxzz"}


def extract() -> dict[str, list[Point]]:
    """Return the six gamma components of open O3, one file each."""
    grid = json_table(DIRECTORY, TABLE_2)
    header = grid[0]
    for column, name in COMPONENTS.items():
        if header[column] != f"\\gamma_{{{name}}}":
            msg = f"unexpected Table 2 header {header[column]!r}"
            raise ValueError(msg)
    if grid[1][0] != "O_3(C_{2v})":
        msg = f"unexpected first block {grid[1][0]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[1:4]:  # CCSD(T) and the two DFT rows of open O3
        for column, name in COMPONENTS.items():
            files.setdefault(f"O3/christodouleas_2009__gamma_{name}.csv", []).append(
                Point(cells[1], "inf", number(cells[column]), "", " | ".join(cells))
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
