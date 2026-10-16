# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Waxler & Cleek (1971): Corning 7940 fused silica, -200 to +20 degrees C.

Table 1 (``tbl_9585f8f840a4``) gives n at ten spectral lines (the header prints each
wavelength in nm with its element, "667.8He", "643.8Cd", ...; the element is split
off) and 23 temperatures, +20 to -200 degrees C in steps of 10. One file per
temperature; rows are the wavelengths in nm as printed. The cells agree with the image
(one misprint, 1.46431 at +10 C and 480.0 nm, is an erratum in the record).

Table 2 (the change of index from +18 to -160 C, with Austin and Pierce's and Rinne's
values) has no quantity and is left to the record's notes.
"""

import re

from _library import Point, number, table, write

DIRECTORY = "Waxler_1971_RefractiveIndicesFused"
TABLE_1 = "tbl_9585f8f840a4"
PREFIX = "SiO2/waxler_1971__"

WAVELENGTHS = ("667.8", "643.8", "587.6", "508.6", "501.6", "480.0", "471.3", "467.8")
WAVELENGTHS += ("435.8", "404.7")
_HEADER = re.compile(r"^(\d+\.\d)(He|Cd|Hg)$")


def _name(temperature: str) -> str:
    """Return a dataset name for a printed temperature ("+20" -> "p20C")."""
    if temperature == "0":
        return "0C"
    sign = {"+": "p", "-": "m"}[temperature[0]]
    return f"{sign}{temperature[1:]}C"


def extract() -> dict[str, list[Point]]:
    """Return the index at each of the 23 temperatures."""
    rows = table(DIRECTORY, TABLE_1)
    matches = [_HEADER.match(cell) for cell in rows[0][1:]]
    wavelengths = tuple(m.group(1) if m else "" for m in matches)
    if wavelengths != WAVELENGTHS or len(rows) != 24:
        msg = f"unexpected Table 1 layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[1:]:
        printed = " | ".join(cells)
        points = [
            Point(wavelength, wavelength, number(cell), "", printed)
            for wavelength, cell in zip(wavelengths, cells[1:], strict=True)
        ]
        files[f"{PREFIX}n_{_name(cells[0])}.csv"] = points
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
