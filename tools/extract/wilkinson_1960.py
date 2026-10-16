# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Wilkinson (1960): the refractivity of nitrogen at 1649-2042 A.

Table II (``tbl_12497a3674f3``; image checked) lists, per emission line, its
intensity, vacuum wavelength (A), emitter, ``(n - 1) x 10^7`` at 0 C and 760 mm, and
the difference observed - calculated (x 10^7) from the paper's Eq. (7). The
refractivities are extracted with the wavelength as row id; the difference column is
the formula's residual and is not extracted (the formula is recorded instead).
Table I (second-order lines, air wavelengths, used for the instrumental shift) holds
no refractivities.
"""

from _library import Point, number, table, write

DIRECTORY = "Wilkinson_1960_RefractiveDispersionNitrogen"
TABLE_II = "tbl_12497a3674f3"


def extract() -> dict[str, list[Point]]:
    """Return Table II's refractivities as a data file."""
    points = []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        wavelength = number(cells[1])
        points.append(
            Point(wavelength, wavelength, number(cells[3]), "", " | ".join(cells))
        )
    return {"N2/wilkinson_1960__experiment_N2.csv": points}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
