# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Bridge & Buckingham (1964): the mean polarizabilities of their Table I.

Table I (``tbl_63a38649ff5e``, CSV) prints, per gas, 10^2 rho_0 (the depolarization
ratio), the mean polarizability alpha (A^3) at 6328 A 'calculated ... from refractive
indices given in the Landolt-Bornstein tables (1962)' and the anisotropy. STRUCTURAL
FIX: the library's cells lost the first column (the gas names) and the column headings;
the image shows the rows in the order Ar, Xe, SF6, H2, N2, O2, CO2, which the script
assigns by position (the seven data rows, image checked value by value). Only the
alpha column is extracted: the depolarization ratios and anisotropies are for a later
phase (and the cells misread several of them, e.g. "0.OO5", "Q.OO1").
"""

from _library import Point, number, table, write

DIRECTORY = "Bridge_1964_PolarizationLaserLight"
TABLE_I = "tbl_63a38649ff5e"
# Row order of the printed table (image); the extracted cells lost the gas column.
GASES = ("Ar", "Xe", "SF6", "H2", "N2", "O2", "CO2")
WAVELENGTH = "6328"  # A, the He-Ne laser line as printed


def extract() -> dict[str, list[Point]]:
    """Return one file per gas with its mean polarizability."""
    rows = table(DIRECTORY, TABLE_I)[1:]
    assert len(rows) == len(GASES), rows
    return {
        f"{gas}/bridge_1964__alpha_{gas}.csv": [
            Point(WAVELENGTH, WAVELENGTH, number(cells[1]), "", " | ".join(cells))
        ]
        for gas, cells in zip(GASES, rows, strict=True)
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
