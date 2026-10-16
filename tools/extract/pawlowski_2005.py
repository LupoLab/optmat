# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""Extract Pawlowski, Jorgensen & Hattig (2005): Cauchy moments of Ne, Ar, Kr.

The tables are read from their JSON cells (Table 2 has no CSV rendering) and were
checked against the images. Columns: k, then S(k) in atomic units from the CCS, CC2,
CCSD and CC3 models, then the "Experiment" columns DOSD (Kumar and Meath) and "Other".
In Tables 3 and 4 the library reads the minus sign of k as "K" ("K2" for -2; the images
print -2): the script restores it. The "Other" cell of k = -2 carries a note mark "c"
(another reference than the rest of the column) and goes to its own file. Table 1 (the
CCSD basis-set convergence for neon) is not extracted.
"""

import re

from _json_table import json_table
from _library import Point, number, write

DIRECTORY = "Pawlowski_2005_CauchyMomentsNe"
TABLES = {"Ne": "tbl_1f9141ba155d", "Ar": "tbl_5e1ed83a0b39", "Kr": "tbl_64847cac1c91"}
COLUMNS = {1: "ccs", 2: "cc2", 3: "ccsd", 4: "cc3", 5: "dosd", 6: "other"}
_ORDER = re.compile(r"^[K-]?(\d+)$")
_MARK = re.compile(r"\^?c$")


def extract() -> dict[str, list[Point]]:
    """Return one file per gas and column (and the k = -2 'Other' value apart)."""
    files: dict[str, list[Point]] = {}
    for gas, table_id in TABLES.items():
        for cells in json_table(DIRECTORY, table_id):
            match = _ORDER.match(cells[0].strip())
            if match is None:  # heading rows
                continue
            order = f"-{match[1]}"
            printed = " | ".join(cells)
            for column, stem in COLUMNS.items():
                cell = cells[column].strip()
                name = "other_s-2" if _MARK.search(cell) else stem
                point = Point(order, order, number(_MARK.sub("", cell)), "", printed)
                files.setdefault(f"{gas}/pawlowski_2005__{name}_{gas}.csv", []).append(
                    point
                )
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
