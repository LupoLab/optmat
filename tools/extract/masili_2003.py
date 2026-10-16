# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Masili & Starace (2003): the dynamic polarizability of helium (Table VI).

Table VI (``tbl_4cd6e269513a``, CSV; image checked, p. 9) prints the dynamic dipole
polarizability alpha(omega) of ground-state helium in atomic units at photon energies
omega = 0.050 to 0.845 a.u., the authors' best coupled-channel result (N_c = N_c' = 6)
in the column "This work", beside other authors' values (not taken). Through the 1s2p
resonance (omega about 0.78) alpha changes sign: the image prints -90.326704 to
-0.336927 at 0.782 to 0.820, and the library read every one of those minus glyphs as a
"2" (290.326704 ... 20.336927; the other columns alike). The rows in ``NEGATIVE`` are
restored here, after checking that each cell is "2" followed by the printed magnitude.
The cell at 0.780, 1765.866727, is positive in the image. The rows below the
resonance (0.050 to 0.770) and from it (0.780 to 0.845) go to two files.
"""

from _library import Point, number, table, write

DIRECTORY = "Masili_2003_StaticDynamicDipole"
# The 1s2p resonance, about 0.7799 a.u.: the rows from 0.780 on are written apart.
RESONANCE = 0.7799
TABLE_VI = "tbl_4cd6e269513a"
# Table VI's omega whose "This work" value the image prints negative, and the magnitude
# it prints (the cell is "2" followed by it).
NEGATIVE = {
    "0.782": "90.326704",
    "0.784": "42.857302",
    "0.785": "33.466465",
    "0.790": "14.990171",
    "0.795": "8.852173",
    "0.800": "5.722310",
    "0.805": "3.768979",
    "0.810": "2.380470",
    "0.815": "1.286429",
    "0.820": "0.336927",
}


def extract() -> dict[str, list[Point]]:
    """Return alpha(omega) of the column "This work", the signs restored."""
    grid = table(DIRECTORY, TABLE_VI)
    if grid[0][:2] != ["v", "This work"]:
        msg = f"unexpected Table VI header {grid[0]!r}"
        raise ValueError(msg)
    points = []
    for cells in grid[1:]:
        omega, cell = cells[0].strip(), cells[1].strip().replace(" ", "")
        value = number(cell)
        if omega in NEGATIVE:
            if cell != "2" + NEGATIVE[omega]:
                msg = f"Table VI at {omega}: {cell!r} is not 2 + the printed magnitude"
                raise ValueError(msg)
            value = "-" + NEGATIVE[omega]
        points.append(Point(omega, omega, value, "", " | ".join(cells)))
    if len(points) != 32 or {p.row for p in points} < set(NEGATIVE):
        msg = f"Table VI: {len(points)} rows"
        raise ValueError(msg)
    below = [p for p in points if float(p.row) < RESONANCE]
    return {
        "He/masili_2003__alpha.csv": below,
        "He/masili_2003__alpha_resonance.csv": points[len(below) :],
    }


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
