# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Burrows et al. (1999): GOME O3 cross sections at selected lines.

Table 1 (``tbl_0aa25c60f740``, CSV; image checked, every cell agrees): room-temperature
cross sections at 253.65 nm (the caption: "in 10^-20 cm^2 molecule^-1") from seven
sources, each with its own temperature (the column headed "T (K)", read ``¹ (K)`` in the
cells; the value column is headed sigma_O3, read ``pO3``). The source names carry their
reference numbers as superscripts, run into the name in the cells ("Hearn12"). Each
source is its own file (a dataset has one temperature).

Table 2 (``tbl_7b7a151c331a``, CSV; image checked, every cell agrees): cross sections at
four visible wavelengths (the column headers, nm) "in 10^-24 cm^2 molecule^-1" from Ref.
26, Ref. 27 and this work, each with its temperature. One file per row.
"""

import re

from _library import Point, table, write

DIRECTORY = "BURROWS_1999_ATMOSPHERICREMOTESENSING"
TABLE_1 = "tbl_0aa25c60f740"
TABLE_2 = "tbl_7b7a151c331a"
PREFIX = "O3/burrows_1999__"
WAVELENGTH_1 = "253.65"  # Table 1's caption

# Table 1's sources (the reference number stripped) and their files.
SOURCES_1 = {
    "Hearn": "hearn_1961",
    "Bass and Paur": "bass_paur_1984",
    "Molina and Molina": "molina_1986",
    "Mauersberger et al": "mauersberger_1987",
    "Yoshino et al": "yoshino_1993",
    "Daumont et al": "daumont_1992",
    "This work": "t293_253",
}
SOURCES_2 = {
    "Ref. 26": "anderson_mauersberger_1992",
    "Ref. 27": "burkholder_talukdar_1994",
    "This work": "t293_visible",
}
_REFERENCE = re.compile(r"\d+$")


def extract() -> dict[str, list[Point]]:
    """Return Table 1 per source and Table 2 per row."""
    files: dict[str, list[Point]] = {}
    for cells in table(DIRECTORY, TABLE_1)[1:]:
        source, temperature, value = (c.strip() for c in cells)
        name = _REFERENCE.sub("", source)  # the reference number run into the name
        float(temperature)
        float(value)
        files[f"{PREFIX}{SOURCES_1[name]}.csv"] = [
            Point(WAVELENGTH_1, WAVELENGTH_1, value, "", " | ".join(cells))
        ]
    header, *rows = table(DIRECTORY, TABLE_2)
    wavelengths = [c.strip() for c in header[2:]]
    for cells in rows:
        source, temperature, *values = (c.strip() for c in cells)
        float(temperature)
        points = []
        for wavelength, value in zip(wavelengths, values, strict=True):
            float(wavelength)
            float(value)
            points.append(Point(wavelength, wavelength, value, "", " | ".join(cells)))
        files[f"{PREFIX}{SOURCES_2[source]}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
