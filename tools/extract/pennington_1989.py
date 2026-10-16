# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Pennington, Henesian & Hellwarth (1989): Table II, argon's c1221 (p. 3007).

Table II gives c1221(-omega, omega, omega, -omega) in 1e-19 esu for the three beam
configurations and the two analyzer angles gamma = 0 and 45 degrees, as ``value \pm
error`` with footnote letters. Only the argon row (configuration (c)) is extracted; the
air rows are not curated. The JSON cells match the table image (checked).
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Pennington_1989_NonlinearIndexAir"
TABLE = "tbl_35e9b208bbb4"
CELL = re.compile(r"^([\d.]+)\\pm([\d.]+)(?:\^\{[a-z,]+\})?$")


def extract() -> dict[str, list[Point]]:
    """Return the argon values, one row per analyzer angle."""
    (argon,) = [cells for cells in json_table(DIRECTORY, TABLE) if cells[1] == "Argon"]
    points = []
    for column, row in ((2, "gamma0"), (3, "gamma45")):
        match = CELL.match(argon[column])
        if match is None:
            raise ValueError(argon[column])
        value, error = match.groups()
        points.append(
            Point(row, "1053", number(value), number(error), " | ".join(argon))
        )
    return {"Ar/pennington_1989__c1221_Ar.csv": points}


def main() -> None:
    """Write the data file."""
    write(extract())


if __name__ == "__main__":
    main()
