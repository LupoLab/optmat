# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Lin et al. (1976): Table I, the Delta alpha of the curated linear molecules.

Table I (p. 819) lists, per molecule, the temperature of the calculation, the
polarizability anisotropy Delta alpha (1e-25 esu, i.e. cm3) with a footnote letter
naming its source (b: Hirschfelder, Curtiss and Bird; c: MacAdam and Ramsey), the
rotational constant B (d, e: Herzberg; not extracted) and the calculated collisionless
n2 (1e-12 esu/mole; not extracted, see the record). The CSV matches the table image
(checked); the footnote letters are printed as trailing characters of the cells.
"""

from _library import Point, number, table, write

DIRECTORY = "Lin_1976_BirefringenceArisingReorientation"
TABLE = "tbl_7ac00c24db21"
GASES = {
    "CO_2": "CO2",
    "N_2O": "N2O",
    "N_2": "N2",
    "O_2": "O2",
    "H_2": "H2",
    "D_2": "D2",
}
SOURCES = {"b": "hirschfelder", "c": "macadam_ramsey"}


def extract() -> dict[str, list[Point]]:
    """Return the reproduced anisotropy of each curated gas (no stated wavelength)."""
    files = {}
    for cells in table(DIRECTORY, TABLE)[1:]:
        if cells[0] not in GASES:
            continue  # CS2, C2H2, HCl, HBr are not curated
        gas = GASES[cells[0]]
        value, letter = cells[2][:-1], cells[2][-1]
        path = f"{gas}/lin_1976__dalpha_{SOURCES[letter]}_{gas}.csv"
        files[path] = [Point("inf", "not stated", number(value), "", " | ".join(cells))]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
