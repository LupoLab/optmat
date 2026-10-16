# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Gill & Heddle (1963): |n - 1| at 1216 A from Rayleigh scattering (Table I).

Table I gives, per gas, the normalised scattering slope S', the depolarisation factor,
the ratio |n0 - 1| / |n0 - 1|_Ar (the square root of their product) and |n0 - 1| x 10^4
(the ratio times argon's value, argon's value being nitrogen's Cerenkov value 5.71e-4
divided by nitrogen's ratio). Both the ratio and the absolute column are extracted,
except xenon's absolute value: |n - 1| = 172e-4 would be a polarizability of 690 a.u.,
outside what the package accepts as a refractivity (and, per the paper and Kingston
1964, not a real refractivity of the atom); xenon is kept as its (excluded) ratio.
"""

from _library import Point, number, table, write

DIRECTORY = "Gill_1963_DeterminationRefractiveIndices"
TABLE = "tbl_8bd54abf02ad"
GASES = {
    "Argon": "Ar",
    "Krypton": "Kr",
    "Xenon": "Xe",
    "Hydrogen": "H2",
    "Nitrogen": "N2",
}
WAVELENGTH = (
    "1216"  # Lyman alpha, as printed in the caption ("Refractivities at 1216 A")
)


def extract() -> dict[str, list[Point]]:
    """Return one ratio file (not for argon, the unit) and one absolute file per gas."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        gas = GASES[cells[0]]
        printed = " | ".join(cells)
        if gas != "Xe":
            files[f"{gas}/gill_1963__absolute_{gas}.csv"] = [
                Point(WAVELENGTH, WAVELENGTH, number(cells[4]), "", printed)
            ]
        if gas != "Ar":
            files[f"{gas}/gill_1963__ratio_{gas}.csv"] = [
                Point(WAVELENGTH, WAVELENGTH, number(cells[3]), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
