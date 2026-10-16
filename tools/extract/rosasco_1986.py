# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Rosasco & Hurst (1986): 2 sigma/24 and 2 B(0)/24 of H2, D2, He (Table 1).

Table 1 (``tbl_1142277dccf1``, published only as JSON cells; the cells agree with the
table image) lists susceptibilities in 10^-19 esu at 1 amagat and 296 K. Only the
electronic 2 sigma/24 (row 1) and the reorientation term 2 B(0)/24 (row 2) are
extracted; the off-resonant Raman sums (rows 4-8) and the D2 Raman amplitude A_0 N_0
(row 9) are in the record's notes. Footnote letters are stripped; a dash is no value.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Rosasco_1986_DispersionElectronicContribution"
TABLE = "tbl_1142277dccf1"
PUMP = "14779"  # cm^-1, Table 1 note a (the probe is 17751 cm^-1)
ROWS = {1: "sigma", 2: "reorientation"}
GASES = {1: "H2", 2: "D2", 3: "He"}
_CELL = re.compile(r"^(?P<value>[\d.]+)\^\{[a-z]\}$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and row."""
    grid = json_table(DIRECTORY, TABLE)
    files = {}
    for row, name in ROWS.items():
        cells = grid[row]
        printed = " | ".join(cells)
        for column, gas in GASES.items():
            if cells[column] == "\N{EM DASH}":
                continue
            match = _CELL.match(cells[column])
            if match is None:
                msg = f"unexpected cell {cells[column]!r}"
                raise ValueError(msg)
            files[f"{gas}/rosasco_1986__{name}_{gas}.csv"] = [
                Point(name, PUMP, number(match["value"]), "", printed)
            ]
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
