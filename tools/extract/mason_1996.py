# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Mason et al. (1996): O3 oscillator strengths of band structures (Table 3).

Table 3 (``tbl_cf055268b405``, cells file only; image checked, every cell agrees):
integrated oscillator strengths of ozone over eight photon-energy intervals, printed
"E_1 (eV) to E_2 (eV)" in the header ("3.675 to 3.956" ...), one row per source:
experiment (Present, Johnstone et al, Davies, Tanaka et al, Brion et al (1984), Inn &
Tanaka, Molina & Molina) and theory (Borowski et al, Tsuneda et al, Hay & Dunning,
Thunemann et al, Nordfors et al, Padial et al). Each cell is printed ``m \times
10^{k}`` with its own power of ten; the value is written in e-notation (``11.2e-5``),
the printed cell kept in ``printed``. "—" is no value.

What the image shows and the cells merge:

- note a on the Present value of the last interval: for that value E_2 = 11.25 eV, not
  the header's 11.30; the row id says so (``10.70-11.25``);
- note b (Brion et al 1984, "Measured at 294 K") and note c (Tsuneda et al: the first
  value CASPT2(g2), the second MRMP) are letters stuck to the cells;
- the Tsuneda cell holds two values one above the other (``8.9 ... 8.0 ...c``): split
  into two rows, ``4.13-6.20 CASPT2(g2)`` and ``4.13-6.20 MRMP``.

The abscissa is the interval's lower edge E_1 in eV (as ``kumar_2010_supplement``), the
row id the interval ``E_1-E_2``.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Mason_1996_VUVOpticalAbsorption"
TABLE_3 = "tbl_cf055268b405"
PREFIX = "O3/mason_1996__"

# The table's source rows and the dataset each goes to.
ROWS = {
    "Present": "present",
    "Johnstone et al": "johnstone_1992",
    "Davies": "davies_1993",
    "Tanaka et al": "tanaka_1953",
    "Brion et al (1984)": "brion_1984",
    "Inn & Tanaka": "inn_tanaka_1953",
    "Molina & Molina": "molina_1986",
    "Borowski et al": "borowski_1995",
    "Tsuneda et al": "tsuneda_1995",
    "Hay & Dunning": "hay_dunning_1977",
    "Thunemann et al": "thunemann_1978",
    "Nordfors et al": "nordfors_1991",
    "Padial et al": "padial_1981",
}
_VALUE = re.compile(r"([\d.]+) \\times 10\^\{(-\d+)\}([abc]?)")
_INTERVAL = re.compile(r"^([\d.]+) to ([\d.]+)$")


def extract() -> dict[str, list[Point]]:
    """Return one file of interval oscillator strengths per source row of Table 3."""
    grid = json_table(DIRECTORY, TABLE_3)
    intervals = []
    for cell in grid[0][1:]:
        match = _INTERVAL.match(cell)
        if match is None:
            msg = f"unexpected interval {cell!r}"
            raise ValueError(msg)
        intervals.append(match.groups())
    files: dict[str, list[Point]] = {}
    for cells in grid[1:]:
        name = cells[0]
        if name in {"Experiment", "Theory"}:
            continue
        points = []
        for (lower, upper), cell in zip(intervals, cells[1:], strict=True):
            if cell == "\N{EM DASH}":
                continue
            values = _VALUE.findall(cell)
            rebuilt = "".join(f"{m} \\times 10^{{{k}}}{n}" for m, k, n in values)
            if not values or rebuilt != cell:
                msg = f"unexpected cell {cell!r}"
                raise ValueError(msg)
            for mantissa, power, note in values:
                edge = "11.25" if note == "a" else upper  # note a: E_2 = 11.25 eV
                row = f"{lower}-{edge}"
                if len(values) == 2:  # Tsuneda: note c, CASPT2(g2) then MRMP
                    row += " CASPT2(g2)" if note == "" else " MRMP"
                printed = f"{name} | {lower} to {upper} | {cell}"
                points.append(
                    Point(
                        row,
                        lower,
                        f"{mantissa}e{power}",
                        "",
                        printed,
                        extra={"x_upper": edge},
                    )
                )
        files[f"{PREFIX}{ROWS[name]}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
