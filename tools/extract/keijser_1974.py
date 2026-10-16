# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Keijser et al. (1974), Physica 76, 585: rotational Raman line broadening.

Table I (``tbl_ef6959254514``, p. 595) exists only as JSON cells, checked against its
image: the broadening coefficients Delta nu_1/2 / rho (HWHM, GHz/amagat) of the S_o(J)
lines at 293 K, "v ± s" cells, "-" for a line not given. The columns are HD (not a
material: the record's notes), normal H2 ("this exp." and "literature") and normal D2.

**Structural fix:** the H2 "literature" cell of S_o(1) holds two values run together,
"0.053 ± 0.001a0.054 ± 0.002b": the image prints them on two lines, marked a (Ref. 17,
Cooper, May and Gupta 1970) and b (Ref. 18, Gupta and May 1972). The script splits the
cell at its note letters and writes one file per reference.

Tables II-VI (populations, cross sections, the theory) are notes in the record.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Keijser_1974_PressureBroadeningRotational"
TABLE_I = "tbl_ef6959254514"
_NOTED = re.compile(r"([\d.]+)\s*±\s*([\d.]+)([ab])")
LITERATURE = {"a": "table1_H2_cooper", "b": "table1_H2_gupta_may"}


def _split(cell: str) -> tuple[str, str]:
    """Return (value, sigma) of a "v ± s" cell."""
    value, _, sigma = cell.partition("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return the H2 and D2 files of Table I."""
    grid = json_table(DIRECTORY, TABLE_I)
    if (
        grid[2][1:] != ["HD", "nH_2", "", "nD_2"]
        or grid[3][2:4] != ["this exp.", "literature"]
        or "GHz/amagat" not in grid[1][1]
    ):
        msg = f"unexpected Table I header {grid[:4]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[4:]:
        label = cells[0].replace("S_o", "S0")  # S_o(J): the lowest vibrational state
        printed = " | ".join(cells)
        for column, name in (
            (2, "H2/keijser_1974__table1_H2.csv"),
            (4, "D2/keijser_1974__table1_D2.csv"),
        ):
            if cells[column] != "-":
                value, sigma = _split(cells[column])
                files.setdefault(name, []).append(
                    Point(label, label, value, sigma, printed)
                )
        for value, sigma, note in _NOTED.findall(cells[3]):
            name = f"H2/keijser_1974__{LITERATURE[note]}.csv"
            files.setdefault(name, []).append(
                Point(label, label, number(value), number(sigma), printed)
            )
        if cells[3] != "-" and not _NOTED.findall(cells[3]):
            msg = f"unexpected literature cell {cells[3]!r}"
            raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
