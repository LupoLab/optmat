# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Chopra & Heddle (1974): Rayleigh scattering of Lyman alpha (Ar, Kr, N2, H2).

- Table 2 (``tbl_6ff9e885d7c3``; image checked, the cells agree): scattering cross
  sections at 1216 A in barns, "Present work" (normalised to Dalgarno et al.'s
  calculated hydrogen value, which is the hydrogen entry), Shardanand and Mikawa
  (1967) and Gill and Heddle (1963). One file per column and gas. The present xenon
  entry, ``1.71 \times 10^{7} \dagger`` (molecular xenon, Xe2), is not a curated
  species and is skipped.
- Table 4 (``tbl_d63bbf88a338``; image checked): 10^4 (n - 1) of argon at 1216 A, one
  file per method row; the Cerenkov value is "at 1210 +- 30 A", which becomes its
  wavelength.

Table 1 (relative cross sections) has no quantity in the format and is described in the
record; Tables 3 (xenon attenuation) and 5 (krypton contributions) are not extracted.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Chopra_1974_PolarizationFreeMeasurements"
TABLE_2 = "tbl_6ff9e885d7c3"
TABLE_4 = "tbl_d63bbf88a338"
GASES = {
    "Hydrogen": "H2",
    "Argon": "Ar",
    "Nitrogen": "N2",
    "Krypton": "Kr",
    "Xenon": "Xe",
}
TABLE_2_COLUMNS = {1: "rayleigh", 2: "shardanand_1967", 3: "gill_1963"}
TABLE_4_ROWS = {
    "Present work": "refractivity",
    "Cerenkov radiation": "heddle_1963",
    "Line displacement": "gashchina_1968",
    "Calculation": "kingston_1964",
}
_AT = re.compile(r"^([\d.]+)\s*\[at\s*(\d+)")


def _table_2() -> dict[str, list[Point]]:
    files = {}
    for cells in table(DIRECTORY, TABLE_2)[1:]:
        gas = GASES[cells[0].strip()]
        for column, stem in TABLE_2_COLUMNS.items():
            cell = cells[column].strip()
            if not cell or "dagger" in cell:
                continue  # blank, or the molecular-xenon entry
            path = f"{gas}/chopra_1974__{stem}_{gas}.csv"
            files[path] = [Point("1216", "1216", number(cell), "", " | ".join(cells))]
    return files


def _table_4() -> dict[str, list[Point]]:
    files = {}
    for cells in table(DIRECTORY, TABLE_4)[1:]:
        label = cells[0].strip().rstrip("\N{DAGGER}\N{DOUBLE DAGGER}\N{SECTION SIGN}")
        cell = cells[1].strip()
        found = _AT.match(cell)
        wavelength, value = (found[2], found[1]) if found else ("1216", cell)
        path = f"Ar/chopra_1974__{TABLE_4_ROWS[label]}_Ar.csv"
        files[path] = [
            Point(wavelength, wavelength, number(value), "", " | ".join(cells))
        ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return the Rayleigh cross sections (Table 2) and the argon refractivities."""
    return _table_2() | _table_4()


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
