# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Peck & Fisher (1964): argon n - 1 at 17 vacuum wavelengths (Tables I, II).

Table I (``tbl_c5791a12e56b``, merged cells, so read from the JSON cells) gives each
measurement reduced to 0 degrees C and to 15 degrees C (with the random error of the
15 degrees C value relative to the 5462 A reference) and the paper's 15 degrees C
formula. Wavelengths above 10 000 A are printed with a thin space (``10 142.573``),
removed here. The reference line (Hg 5462.258 A) prints ``...`` for its relative
error: it is the reference of the relative measurements, so its relative error is 0
by construction (its absolute uncertainty is the dataset's systematic one).

Table II compares the paper's formula with Larsen's (1934) measurements at 0 degrees C;
Larsen's column is a reproduction.
"""

import json

from _library import LIBRARY, Point, number, table, write

DIRECTORY = "Peck_1964_DispersionArgon"
TABLE_I = "tbl_c5791a12e56b"
TABLE_II = "tbl_1681f4cff0cd"


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def extract() -> dict[str, list[Point]]:
    """Return Table I's three columns and Table II's two as data files."""
    measured_15, measured_0, fit_15 = [], [], []
    for cells in grid(TABLE_I)[1:18]:  # the last row is the rms-deviation note
        printed = " | ".join(cells)
        wavelength = cells[1].replace(" ", "")
        value_15, _, error = cells[3].replace("...", "").partition("±")
        sigma = "0" if error == "" else number(error.rstrip("a"))
        measured_15.append(
            Point(wavelength, wavelength, number(value_15), sigma, printed)
        )
        measured_0.append(Point(wavelength, wavelength, number(cells[2]), "", printed))
        fit_15.append(Point(wavelength, wavelength, number(cells[4]), "", printed))
    fit_0, larsen = [], []
    for cells in table(DIRECTORY, TABLE_II)[1:]:
        printed = " | ".join(cells)
        fit_0.append(Point(cells[0], cells[0], number(cells[1]), "", printed))
        larsen.append(Point(cells[0], cells[0], number(cells[2]), "", printed))
    return {
        "Ar/peck_1964__measured_15C.csv": measured_15,
        "Ar/peck_1964__measured_0C.csv": measured_0,
        "Ar/peck_1964__fit_15C.csv": fit_15,
        "Ar/peck_1964__fit_0C_table2.csv": fit_0,
        "Ar/peck_1964__larsen_1934.csv": larsen,
    }


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
