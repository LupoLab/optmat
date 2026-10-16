# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Rodney & Spindler (1954): fused-quartz glass, seven prisms, 24 degrees C.

The library copy is a scan with an OCR text layer; every table was checked against the
page images (pdf pages 3 and 4, journal pp. 187-188).

Structural fix (captions): the extraction attaches Table 2's caption ("Most probable
values of refractive indices of pure fused-quartz glass at 24 C ... Computed from eq
(1)") to the residual table ``tbl_947a659cc842`` and turns Table 1's caption into a
table note, while the grid ``tbl_39b3d33089ff`` has no caption. The page image shows
``tbl_947a659cc842`` under the caption of Table 1 ("Computed indices of refraction of
fused-quartz glass (Heraeus B) and residuals (o-c) x 10^5 for temperature of 24 C")
and ``tbl_39b3d33089ff`` under that of Table 2. The script uses the image's
assignment.

- **Table 1** (``tbl_947a659cc842``, cells only): 24 wavelengths (um; the leading
  zero is omitted in print, ".36117"), the index computed from Eq. (1), and the
  residuals (o - c) x 10^5 of seven prisms against Eq. (1): Heraeus A, B, C, Nieder
  disk, General Electric G1, G2, Corning silica glass. A dashed cell ("-" in the
  cells) is a line not measured on that prism and is skipped. The cells agree with the
  image (the header "G?-c" of the image is G2-c).
- **Table 2** (``tbl_39b3d33089ff``): n from Eq. (1) at 0.34-3.50 um, printed in two
  column pairs (0.34-0.66, 0.67-3.50); rows in that order.
- **Table 3** (``tbl_69b8ae37c6ed``, cells only): 10^7 x Delta n / Delta t and n at 15,
  20, 25, 30 and 35 C, 0.40-0.70 um. The 30 C value at 0.70 um is printed "0.45541"
  (image checked): a misprint, corrected by an erratum in the record, not here.

The x column and the row id are the wavelength with the omitted leading zero restored
(".36117" -> "0.36117"); the printed column keeps the cells as printed.
"""

import re

from _json_table import json_table
from _library import Point, number, table, write

DIRECTORY = "Rodney_1954_IndexRefractionFused"
TABLE_1 = "tbl_947a659cc842"
TABLE_2 = "tbl_39b3d33089ff"
TABLE_3 = "tbl_69b8ae37c6ed"
PREFIX = "SiO2/rodney_1954__"

# Table 1 column -> dataset id.
TABLE_1_COLUMNS = {
    2: "computed_table1",
    3: "residual_heraeus_a",
    4: "residual_heraeus_b",
    5: "residual_heraeus_c",
    6: "residual_nieder",
    7: "residual_ge_g1",
    8: "residual_ge_g2",
    9: "residual_corning",
}
TABLE_1_HEADER = ("A-c", "B-c", "C-c", "", "G1-c", "G2-c", "")
TABLE_3_TEMPERATURES = ("15", "20", "25", "30", "35")
_DEGREES = re.compile(r"t = \{(\d+)\}")


def _wavelength(cell: str) -> str:
    """Return a printed wavelength (um) with its omitted leading zero restored."""
    text = cell.strip()
    return number("0" + text if text.startswith(".") else text)


def _table_1() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_1)
    if tuple(rows[1][3:]) != TABLE_1_HEADER or rows[2][1] != "μ":
        msg = f"unexpected Table 1 header {rows[:3]!r}"
        raise ValueError(msg)
    body = rows[3:]
    if len(body) != 24:
        msg = f"Table 1 has {len(body)} wavelengths, not 24"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in body:
        wavelength = _wavelength(cells[1])
        printed = " | ".join(cells)
        for column, dataset in TABLE_1_COLUMNS.items():
            if cells[column] == "-":  # a dashed cell: not measured on this prism
                continue
            point = Point(wavelength, wavelength, number(cells[column]), "", printed)
            files.setdefault(f"{PREFIX}{dataset}.csv", []).append(point)
    return files


def _table_2() -> dict[str, list[Point]]:
    rows = table(DIRECTORY, TABLE_2)
    if rows[0] != ["Wave-length", "Index", "Wave-length", "Index"] or rows[1][0] != "μ":
        msg = f"unexpected Table 2 header {rows[:2]!r}"
        raise ValueError(msg)
    points: list[Point] = []
    for half in (0, 2):
        for cells in rows[2:]:
            if not cells[half]:
                continue
            wavelength = _wavelength(cells[half])
            printed = f"{cells[half]} | {cells[half + 1]}"
            value = number(cells[half + 1])
            points.append(Point(wavelength, wavelength, value, "", printed))
    if len(points) != 65:
        msg = f"Table 2 has {len(points)} wavelengths, not 65"
        raise ValueError(msg)
    return {f"{PREFIX}table2.csv": points}


def _table_3() -> dict[str, list[Point]]:
    rows = json_table(DIRECTORY, TABLE_3)
    degrees = tuple(_DEGREES.match(cell).group(1) for cell in rows[0][2:])  # pyright: ignore[reportOptionalMemberAccess]
    if degrees != TABLE_3_TEMPERATURES or len(rows) != 33:
        msg = f"unexpected Table 3 layout {rows[0]!r}, {len(rows)} rows"
        raise ValueError(msg)
    files: dict[str, list[Point]] = {}
    for cells in rows[2:]:
        wavelength = _wavelength(cells[0])
        printed = " | ".join(cells)
        point = Point(wavelength, wavelength, number(cells[1]), "", printed)
        files.setdefault(f"{PREFIX}dndt.csv", []).append(point)
        for temperature, cell in zip(TABLE_3_TEMPERATURES, cells[2:], strict=True):
            point = Point(wavelength, wavelength, number(cell), "", printed)
            files.setdefault(f"{PREFIX}table3_{temperature}C.csv", []).append(point)
    return files


def extract() -> dict[str, list[Point]]:
    """Return Table 1's computed index and residuals, Table 2, and Table 3."""
    return {**_table_1(), **_table_2(), **_table_3()}


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
