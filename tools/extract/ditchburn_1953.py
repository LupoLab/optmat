# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Ditchburn & Heddle (1953): O2 continuum absorption coefficients, 1360-1670 A.

Table 1 (``tbl_9b71613fa0c5``, cells file only; image checked, the cells agree cell for
cell) lists for about fifty hydrogen emission lines the wave number nu (cm^-1, column 1)
and the absorption coefficient alpha (cm^-1 at 1 atm, 0 C) measured with the 205 cm path
(column 2) and the 0.043 cm path (column 3), each "<alpha> (<number of observations>)"
or a dash, then, on the last line of each group of neighbouring lines, the group's mean
wave number (column 4), mean alpha "<value> +- <error>" (column 5) and number of
observations (column 6). Two header rows (the column numbers and names, then the two
path lengths under the brace over columns 2-3). The image separates the groups by extra
space; each mean stands on the group's last line. The script checks that the
observations of each group's lines add up to column 6 (they do, for all 21 groups).

Files: the group means with their errors (``means``, sigma from column 5), and the
single-line values of the long path (``long_path``) and the short path (``short_path``);
the number of observations stays in the printed column. Row ids are the printed wave
numbers (unique within each file; 66780 appears once per path).
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Ditchburn_1953_AbsorptionCrossSections"
TABLE_1 = "tbl_9b71613fa0c5"
PREFIX = "O2/ditchburn_1953__"

_LINE = re.compile(r"^(\d+) \((\d+)\)$")
_MEAN = re.compile(r"^(\d+) ± (\d+)$")
_DASH = "\N{EM DASH}"


def _line(cell: str) -> tuple[str, int] | None:
    """Return a single-line value and its number of observations, or None for a dash."""
    if cell == _DASH:
        return None
    match = _LINE.match(cell)
    if match is None:
        msg = f"unexpected alpha cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), int(match.group(2))


def extract() -> dict[str, list[Point]]:
    """Return the group means and the single-line values of the two paths."""
    rows = json_table(DIRECTORY, TABLE_1)
    if rows[1][1:3] != ["l = 205", "l = 0\\cdot043"]:
        msg = f"unexpected path header {rows[1]!r}"
        raise ValueError(msg)
    means, long_path, short_path = [], [], []
    observations = 0
    for cells in rows[2:]:
        nu, long_cell, short_cell, nu_mean, alpha_mean, count = (
            c.strip() for c in cells
        )
        int(nu)
        printed = " | ".join((nu, long_cell, short_cell, nu_mean, alpha_mean, count))
        for cell, points in ((long_cell, long_path), (short_cell, short_path)):
            value = _line(cell)
            if value is not None:
                points.append(Point(nu, nu, value[0], "", printed))
                observations += value[1]
        if nu_mean:
            match = _MEAN.match(alpha_mean)
            if match is None:
                msg = f"unexpected mean cell {alpha_mean!r}"
                raise ValueError(msg)
            if observations != int(count):
                msg = f"group ending {nu}: {observations} observations, not {count}"
                raise ValueError(msg)
            observations = 0
            int(nu_mean)
            means.append(
                Point(nu_mean, nu_mean, match.group(1), match.group(2), printed)
            )
    if observations:
        msg = "lines after the last group mean"
        raise ValueError(msg)
    return {
        f"{PREFIX}means.csv": means,
        f"{PREFIX}long_path.csv": long_path,
        f"{PREFIX}short_path.csv": short_path,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
