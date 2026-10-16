# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Heddle, Jennings & Parsons (1963): Ar and N2 n - 1 by the Cerenkov method.

Table I (``tbl_2c6e9aed4bea``, merged header cells, so read from the JSON cells) has a
nitrogen part (rows 3-5) and an argon part (rows 9-11): the effective wavelength with
the detector band after ``+-``, the measured refractivity (x 10^5, STP) with its
standard deviation, and the value from Wilkinson's (N2) or Larsen's (Ar) dispersion
formula (``...`` where none).
"""

import json

from _library import LIBRARY, Point, number, write

DIRECTORY = "Heddle_1963_DeterminationRefractiveIndices"
TABLE = "tbl_2c6e9aed4bea"
PARTS = {"N2": (slice(3, 6), "wilkinson"), "Ar": (slice(9, 12), "larsen")}


def grid(table_id: str) -> list[list[str]]:
    """Return a library table's cell grid from its JSON file (merged cells once)."""
    path = LIBRARY / DIRECTORY / "tables" / f"{table_id}.json"
    cells = json.loads(path.read_text(encoding="utf-8"))
    rows = [[""] * cells["columns"] for _ in range(cells["rows"])]
    for cell in cells["cells"]:
        rows[cell["row"]][cell["column"]] = cell["text"].strip()
    return rows


def extract() -> dict[str, list[Point]]:
    """Return the measured and the formula column of each gas."""
    rows = grid(TABLE)
    files = {}
    for gas, (part, formula) in PARTS.items():
        measured, calculated = [], []
        for cells in rows[part]:
            printed = " | ".join(cells)
            wavelength = cells[1].split("±")[0]
            value, _, sigma = cells[2].partition("±")
            measured.append(
                Point(wavelength, wavelength, number(value), sigma, printed)
            )
            if cells[3] != "...":
                calculated.append(
                    Point(wavelength, wavelength, number(cells[3]), "", printed)
                )
        files[f"{gas}/heddle_1963__experiment_{gas}.csv"] = measured
        files[f"{gas}/heddle_1963__{formula}_formula_{gas}.csv"] = calculated
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
