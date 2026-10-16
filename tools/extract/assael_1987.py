# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Assael, Mixafendi & Wakeham (1987): the zero-density viscosity of normal D2.

Table 2 (``tbl_8b226138c853``, p. 191) prints eta/(uPa s) against T/K from 30 K to
2200 K in four column pairs, read here column by column into one file, sorted by
temperature. Four cells of the extracted table differ from the page, checked against a
400 dpi rendering of the page image (2026-09-30), each also fitting its neighbours:

- the first temperature of the third pair is 1080, not 1060 (1080, 1090, 1100 follow);
- 280 K is 12.06, not 12.08 (steps of 0.29 from 11.77 and to 12.35);
- 840 K is 25.68, not 25.69 (25.46, 25.68, 25.89);
- 1200 K is 32.99, not 32.90 (32.80, 32.99, 33.19).

Every other cell agrees with the image. The page prints 1300 K as 34.90, one unit in the
last digit off the smooth run (34.72, 34.90, 35.10); it is kept as printed.
"""

from _library import Point, number, table, write

DIRECTORY = "Assael_1987_ViscosityNormalDeuterium"
TABLE_2 = "tbl_8b226138c853"
FILE = "D2/assael_1987__table2.csv"
FIXES = {
    # (row index in the extracted table, column): (extracted, printed on the page)
    (1, 4): ("1060", "1080"),
    (33, 1): ("12.08", "12.06"),
    (33, 3): ("25.69", "25.68"),
    (13, 5): ("32.90", "32.99"),
}


def extract() -> dict[str, list[Point]]:
    """Return Table 2's viscosities, one point per temperature."""
    grid = table(DIRECTORY, TABLE_2)
    if grid[0][0] != "T/K" or "Pa s" not in grid[0][1]:
        msg = f"unexpected Table 2 header {grid[0]!r}"
        raise ValueError(msg)
    for (row, column), (extracted, printed) in FIXES.items():
        if grid[row][column] != extracted:
            msg = f"Table 2 cell {(row, column)} is {grid[row][column]!r}"
            raise ValueError(msg)
        grid[row][column] = printed
    points = []
    for cells in grid[1:]:
        for i in range(0, len(cells), 2):
            if cells[i]:
                temperature, value = number(cells[i]), number(cells[i + 1])
                points.append(
                    Point(
                        temperature, temperature, value, "", f"{temperature} | {value}"
                    )
                )
    points.sort(key=lambda p: float(p.x))
    return {FILE: points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
