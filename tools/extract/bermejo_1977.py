# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bermejo, Escribano and Orza (1977): absolute Raman intensities of methane.

Table I (``tbl_602282e66a5e``, p. 348, continued in ``tbl_db3b249551b4``, p. 349; image
checked) prints, per band, the band position (``a``: 'Our measurement'), the absolute
intensity I = g(45 a'^2 + 7 g'^2) in units of N x 10^-32 cm^4 g^-1 (N Avogadro's
number: angstrom^4/amu) and the depolarization ratio. Only CH4 is curated: its block
is the left half of the first part's rows 2-6. The deuterated molecules (CH3D, CH2D2,
CHD3, CD4) are isotopologues, left to the record's notes.

From CH4's block, the band positions (``raman_shift``) and the intensities
(``raman_scattering_activity``) are taken. Not written: the nu4 intensity
``\leqslant 0.05^{Q}`` (an upper limit) and every CH4 depolarization ratio (``0^{d}``,
``3/4^{d}``: 'Expected values, not measured'). The nu3 position is printed as two
braced maxima, ``{3016.7^a 3021.3^a``; the cells join them in one cell, split here.
The 2nu2 position is printed ``~3071``: the value 3071 is written and the ``~`` kept in
``printed``.

Table II (``tbl_29036d066e28``, p. 350; image checked): the calculated (from the fitted
alpha^S tensors) and the adopted experimental intensities, the latter as
``value±sigma``. CH4 is rows 1-4; nu4's ``≤0.24`` (an upper limit) is not written.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Bermejo_1977_AbsoluteRamanIntensities"
TABLE_I = "tbl_602282e66a5e"
TABLE_II = "tbl_29036d066e28"
# Table I's CH4 rows (left half) -> row label; nu4 has only an upper limit.
BANDS = {
    2: "nu1 (a1)",
    3: "2nu2 (A1+E)",
    4: "nu2 (e)",
    5: "nu3 (f2)",
    6: "nu4 (f2)",
}
TABLE_II_BANDS = {1: "nu1", 2: "nu2", 3: "nu3"}
_MARK = re.compile(r"\^\{[a-zQ]\}")


def _strip(cell: str) -> str:
    """Return a cell without its footnote marks and the approximate sign."""
    return _MARK.sub("", cell).replace(r"\sim", "").strip()


def extract() -> dict[str, list[Point]]:
    """Return CH4's rows of Tables I and II as data files."""
    grid = json_table(DIRECTORY, TABLE_I)
    if grid[1][0] != "CH_{4}":
        msg = f"Table I row 1 is {grid[1][0]!r}, not CH4"
        raise ValueError(msg)
    shifts: list[Point] = []
    activities: list[Point] = []
    for row, band in BANDS.items():
        cells = grid[row][:4]
        printed = " | ".join(cells)
        position = cells[1]
        if position.startswith("\\{"):
            # The braced pair of maxima of nu3, joined in one cell: split at the marks.
            first, second = _strip(
                position.removeprefix("\\{").replace("^{a}", " ", 1)
            ).split()
            shifts.append(
                Point(
                    f"{band} maximum 1", f"{band} maximum 1", number(first), "", printed
                )
            )
            shifts.append(
                Point(
                    f"{band} maximum 2",
                    f"{band} maximum 2",
                    number(second),
                    "",
                    printed,
                )
            )
        else:
            shifts.append(Point(band, band, number(_strip(position)), "", printed))
        if row != 6:  # nu4: an upper limit
            activities.append(Point(band, band, number(_strip(cells[2])), "", printed))
    grid2 = json_table(DIRECTORY, TABLE_II)
    if grid2[1][0] != "CH_4":
        msg = f"Table II row 1 is {grid2[1][0]!r}, not CH4"
        raise ValueError(msg)
    calculated: list[Point] = []
    experimental: list[Point] = []
    for row, band in TABLE_II_BANDS.items():
        cells = grid2[row]
        printed = " | ".join(cells)
        calculated.append(Point(band, band, number(cells[2]), "", printed))
        value, sigma = cells[3].split("±")
        experimental.append(Point(band, band, number(value), number(sigma), printed))
    return {
        "CH4/bermejo_1977__shift_CH4.csv": shifts,
        "CH4/bermejo_1977__intensity_CH4.csv": activities,
        "CH4/bermejo_1977__intensity_exp_CH4.csv": experimental,
        "CH4/bermejo_1977__intensity_calc_CH4.csv": calculated,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
