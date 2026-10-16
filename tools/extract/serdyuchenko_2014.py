# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Serdyuchenko et al. (2014): O3 cross sections near 220 K, Chappuis and Wulf.

Table 4 (``tbl_75e505961474``, CSV; image checked, every cell agrees): cross sections
("x 10^-22 cm^2 molecule^-1", the caption) at 13 wavelengths (543.667-1046.766 nm,
"taken from El Helou et al. (2005)"), from El Helou et al. (2005) at 222-223 K and at
150 K, Burkholder and Talukdar (1994) at 220 K, Bogumil et al. (2003) at 223 K, BMD at
218 K and this work at 223 K. The two header rows are split across cells in the CSV
("El Helou et a" | "l. (2005)"); the image shows one heading "El Helou et al. (2005)"
over the two temperature columns. "-" is no value. One file per column.
"""

from _library import Point, table, write

DIRECTORY = "Serdyuchenko_2014_HighSpectralResolution"
TABLE_4 = "tbl_75e505961474"
PREFIX = "O3/serdyuchenko_2014__"
COLUMNS = (
    "el_helou_2005_223",
    "el_helou_2005_150",
    "burkholder_talukdar_1994_220",
    "bogumil_2003_223",
    "bmd_218",
    "this_work_223",
)


def extract() -> dict[str, list[Point]]:
    """Return one file per source column of Table 4."""
    files: dict[str, list[Point]] = {f"{PREFIX}{c}.csv": [] for c in COLUMNS}
    for raw in table(DIRECTORY, TABLE_4)[2:]:
        cells = [c.strip() for c in raw]
        wavelength = cells[0]
        float(wavelength)
        printed = " | ".join(cells)
        for name, cell in zip(COLUMNS, cells[1:], strict=True):
            if cell == "-":
                continue
            float(cell)
            files[f"{PREFIX}{name}.csv"].append(
                Point(wavelength, wavelength, cell, "", printed)
            )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
