# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rusch (1923): argon refractivity in the visible and at 6.76 um.

The second Tabelle 7 (``tbl_1f78cf25d43e``, cells file; image checked) lists for argon,
per wavelength ``lambda x 10^7 cm`` (i.e. in nm), ``(nu_0 - 1) x 10^7``: observed
("beob."), calculated with the constants fitted to the visible values only ("ber. I")
and to all five ("ber. II"), and the two differences (not extracted). German decimal
commas are read as points; the printed cells are kept. Tabelle 5 repeats the observed
and "ber. II" columns (and prints 576,95 as "567,95", image checked); it is not
extracted again.

The first row, 6763,2, is the infrared value (the centroid of the calcite reststrahlen,
a vacuum wavelength); the other four are the Hg and H lines (printed air wavelengths).
They go to separate files because their wavelength media differ.

Air (Tabellen 3, 4, 7 first part, 8, 9) is not a curated species; Tabellen 1 and 2 are
calcite reststrahlen wavelengths, and Tabelle 6 holds dispersion constants (recorded as
formulas in the source record).
"""

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rusch_1923_UberLichtbrechungGasen"
TABLE_7_ARGON = "tbl_1f78cf25d43e"
INFRARED = "6763,2"
COLUMNS = {1: "experiment", 2: "fit_I", 3: "fit_II"}


def _german(cell: str) -> str:
    """Return a printed German number (decimal comma) as a plain number string."""
    return number(cell.strip().replace(",", "."))


def extract() -> dict[str, list[Point]]:
    """Return the observed and calculated argon columns, visible and infrared."""
    files: dict[str, list[Point]] = {}
    for cells in json_table(DIRECTORY, TABLE_7_ARGON)[2:]:
        wavelength = _german(cells[0])
        infrared = cells[0].strip() == INFRARED
        printed = " | ".join(cells)
        for column, stem in COLUMNS.items():
            name = f"{stem}_ir" if infrared else stem
            path = f"Ar/rusch_1923__{name}_Ar.csv"
            point = Point(wavelength, wavelength, _german(cells[column]), "", printed)
            files.setdefault(path, []).append(point)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
