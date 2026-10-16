# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Russell & Roh (1987), J. Mol. Spectrosc. 124, 240: D2 Q-branch widths.

Table I (``tbl_99029db61edc``, JSON cells only; checked against the table image, p.
241): the pressure-broadening coefficient a (10^-3 cm-1/amagat, the FWHM law of Eq.
(1)) of D2 Q(0)-Q(4) from this study, and the values of Murray & Javan and of
Robertson & Hunt printed beside them ("-" where they give none). The diffusion
coefficients D_0 and the rms errors delta are in the record (formula and notes), not
written here. Cells print "v ± s".
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Russell_1987_HighResolutionCARS"
TABLE_I = "tbl_99029db61edc"
COLUMNS = {1: "table1_a", 4: "table1_murray_javan_a", 7: "table1_robertson_hunt_a"}


def _split(cell: str) -> tuple[str, str]:
    """Return (value, sigma) of a "v ± s" cell."""
    value, _, sigma = cell.replace(" ", "").partition("±")
    return number(value), number(sigma)


def extract() -> dict[str, list[Point]]:
    """Return the a columns of Table I."""
    grid = json_table(DIRECTORY, TABLE_I)
    heads = [
        [
            "Transition",
            "This Study",
            "",
            "\\delta",
            "Murray & Javan",
            "",
            "\\delta",
            "Robertson & Hunt",
        ],
        ["", "a", "D_{o}", "", "a", "D_{o}", "", "a"],
    ]
    if grid[:2] != heads:
        msg = f"unexpected Table I heads {grid[:2]!r}"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in grid[2:]:
        label = cells[0]
        for column, name in COLUMNS.items():
            if cells[column] == "-":  # no value from that source
                continue
            value, sigma = _split(cells[column])
            files.setdefault(f"D2/russell_1987__{name}.csv", []).append(
                Point(label, label, value, sigma, " | ".join(cells))
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
