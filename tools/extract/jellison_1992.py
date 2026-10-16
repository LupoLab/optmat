# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Jellison (1992): n and k of (100) silicon (Table 1) and of the other faces.

Table 1 (two pages, ``TABLE_1``: pp. 3-4 of the pdf) prints, for (100) silicon at
room temperature, 187 rows of wavelength (nm), photon energy (eV), n and k, each with
its error limit (``1.688±0.015``), in two side-by-side column blocks per page: 234-338
and 340-444 nm (p. 3), 446-640 and 645-840 nm (p. 4). The wavelength is the abscissa
(the measurement grid: 2 nm steps to 450 nm, 5 nm above); the energy is kept in
``printed``. The value and its error limit go to ``value`` and ``sigma``.

Table 2 (``TABLE_2``, p. 6; no CSV, read from its JSON cells; image checked) compares
at seven wavelengths the (100), (111) and (110) faces of this work (ref. a) with Palik's
Handbook (b) and Jellison & Modine 1982 (c). The (111) and (110) rows of this work are
extracted (checks); its (100) rows are Table 1's, which the script verifies; the b and c
rows are other papers' values and are not extracted.

Both pages of Table 1 were compared with their images. The image prints three
misprints, which the library reads as printed and the record corrects by errata: 595.0
nm n "2.948" (3.948), 705.0 nm n "2.766" (3.766), 495.0 nm error limit of n "0.002"
(0.020). The energy printed at 386.0 nm, "3.221", is 3.211 (the column is 1239.5/lambda
in every other row); it is not extracted.
"""

import re

from _json_table import json_table
from _library import Point, table, write

DIRECTORY = "Jellison_1992_OpticalFunctionsSilicon"
TABLE_1 = ("tbl_290c9b5e3da6", "tbl_6b86dd05f9d8")
TABLE_2 = "tbl_364290431ef8"
HEADER = ["Wavelength (nm)", "Energy (eV)", "n", "k"] * 2
ROWS = 187
PREFIX = "Si/jellison_1992__"

_WITH_ERROR = re.compile(r"^(\d+\.\d+)±(\d+\.\d+)$")


def _split(cell: str) -> tuple[str, str]:
    """Return the value and the error limit of a cell ``value±error``."""
    match = _WITH_ERROR.match(cell.strip())
    if not match:
        msg = f"unexpected cell {cell!r}"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def _table_1() -> list[list[str]]:
    """Return Table 1's rows (wavelength, energy, n, k) in wavelength order."""
    rows: list[list[str]] = []
    for table_id in TABLE_1:
        cells = table(DIRECTORY, table_id)
        if cells[0] != HEADER:
            msg = f"{table_id}: unexpected header {cells[0]}"
            raise ValueError(msg)
        for row in cells[1:]:
            rows += [half for half in (row[:4], row[4:]) if half[0].strip()]
    rows = [[cell.strip() for cell in row] for row in rows]
    rows.sort(key=lambda row: float(row[0]))
    if len(rows) != ROWS or len({row[0] for row in rows}) != ROWS:
        msg = f"Table 1: {len(rows)} rows, not {ROWS} distinct wavelengths"
        raise ValueError(msg)
    return rows


def _points(rows: list[list[str]], column: int) -> list[Point]:
    points = []
    for row in rows:
        value, error = _split(row[column])
        points.append(Point(row[0], row[0], value, error, " | ".join(row)))
    return points


def _table_2(table_1: list[list[str]]) -> dict[str, list[Point]]:
    """Return the (111) and (110) faces of this work from Table 2."""
    grid = json_table(DIRECTORY, TABLE_2)
    if grid[0] != ["Wavelength (nm)", "Energy (eV)", "face", "n", "k", "Ref."]:
        msg = f"Table 2: unexpected header {grid[0]}"
        raise ValueError(msg)
    by_wavelength = {f"{float(row[0]):.0f}": row for row in table_1}
    faces: dict[str, list[list[str]]] = {"111": [], "110": []}
    wavelength = energy = ""
    for row in grid[1:]:
        wavelength, energy = row[0] or wavelength, row[1] or energy
        face, n, k, ref = row[2:]
        if ref != "a":
            continue
        cells = [wavelength, energy, face, n, k, ref]
        if face == "100":
            # Table 1's row, to its energy and every digit.
            if by_wavelength[wavelength][1:] != [energy, n, k]:
                msg = f"Table 2 {wavelength} nm (100) differs from Table 1"
                raise ValueError(msg)
            continue
        faces[face].append(cells)
    files = {}
    for face, rows in faces.items():
        if len(rows) != 7:
            msg = f"Table 2: {len(rows)} rows of the ({face}) face, not 7"
            raise ValueError(msg)
        for name, column in (("n", 3), ("k", 4)):
            files[f"{PREFIX}{name}_{face}.csv"] = [
                Point(r[0], r[0], *_split(r[column]), " | ".join(r)) for r in rows
            ]
    return files


def extract() -> dict[str, list[Point]]:
    """Return n and k of (100) silicon, and of the (111) and (110) faces."""
    rows = _table_1()
    return {
        f"{PREFIX}n_100.csv": _points(rows, 2),
        f"{PREFIX}k_100.csv": _points(rows, 3),
    } | _table_2(rows)


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
