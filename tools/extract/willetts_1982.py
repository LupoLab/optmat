# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

r"""Extract Willetts and Harris (1982): He, N2 and CO2 refractivity at 10.59 um.

Table 2 (``tbl_5618cf9a9efd``, JSON only; image checked) gives (n - 1) x 10^6 and the
Gladstone-Dale constant for He, N2, CO2 and He:N2:CO2 mixtures at 10.59 um. The three
pure gases' measured refractivities are extracted; the mixtures (a check of additivity,
three of them calculated) and K are not.
"""

import re

from _json_table import json_table
from _library import Point, write

DIRECTORY = "Willetts_1982_InvestigationIntoOrigin"
TABLE = "tbl_5618cf9a9efd"
# The volumetric ratio He:N2:CO2 of the pure gases' rows.
PURE = {"1 : 0 : 0": "He", "0 : 1 : 0": "N2", "0 : 0 : 1": "CO2"}
_VALUE = re.compile(r"^\s*([\d.]+)\s*±\s*([\d.]+)\s*$")


def extract() -> dict[str, list[Point]]:
    """Return the three pure gases' refractivities as data files."""
    files: dict[str, list[Point]] = {}
    rows = json_table(DIRECTORY, TABLE)
    measured = rows.index(next(r for r in rows if r[0].startswith("(a)")))
    calculated = rows.index(next(r for r in rows if r[0].startswith("(b)")))
    for cells in rows[measured + 1 : calculated]:
        gas = PURE.get(cells[0])
        if gas is None:
            continue  # the 4:2:1 mixture
        match = _VALUE.match(cells[1])
        if match is None:
            msg = f"no value in {cells[1]!r}"
            raise ValueError(msg)
        files[f"{gas}/willetts_1982__measured_{gas}.csv"] = [
            Point("10.59", "10.59", match.group(1), match.group(2), " | ".join(cells))
        ]
    if len(files) != len(PURE):
        msg = f"expected {len(PURE)} pure gases, found {sorted(files)}"
        raise ValueError(msg)
    return files


def main() -> None:
    """Write the data files."""
    write(extract())


if __name__ == "__main__":
    main()
